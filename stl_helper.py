import numpy as np
import struct
from collections import defaultdict

class STLHelper:
    def find_openings(self, stl_path):
        """
        Advanced Hole Detector:
        1. Checks for naked edges (Surface Meshes).
        2. If none, checks for 'Planar Vents' on boundary faces (Watertight Solids).
        """
        triangles = self._load_stl(stl_path)
        if not triangles: return []
        np_tris = np.array(triangles)

        # --- METHOD 1: Naked Edges (For Surface/Shell meshes) ---
        naked_centers = self._find_naked_edge_holes(triangles)
        if naked_centers:
            return naked_centers

        # --- METHOD 2: Internal Planar Loops (For solid Watertight walls with vents) ---
        planar_centers = self._find_planar_holes(np_tris)
        return planar_centers

    def _find_naked_edge_holes(self, triangles):
        edges = defaultdict(int)
        for tri in triangles:
            # Rounding to avoid floating point mismatch
            t = [tuple(np.round(p, 4)) for p in tri]
            pairs = [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])]
            for p1, p2 in pairs:
                edge = tuple(sorted((p1, p2)))
                edges[edge] += 1

        naked_edges = [e for e, count in edges.items() if count == 1]
        if not naked_edges: return []

        # Cluster the edges
        points = set()
        for e in naked_edges:
            points.add(e[0])
            points.add(e[1])

        pts = np.array(list(points))

        # Fallback clustering
        try:
            from sklearn.cluster import DBSCAN
            clustering = DBSCAN(eps=0.5, min_samples=3).fit(pts)
            centers = []
            for label in set(clustering.labels_):
                if label == -1: continue
                cluster_pts = pts[clustering.labels_ == label]
                center = np.mean(cluster_pts, axis=0)
                
                # Approximate area for naked edges
                mins = np.min(cluster_pts, axis=0)
                maxs = np.max(cluster_pts, axis=0)
                dims = sorted(maxs - mins)
                area = dims[1] * dims[2]
                if area < 0.01: area = 0.25 # Fallback 0.25 m2 if calculation is too small
                
                # FIXED: Pack area_m2 into dictionary
                centers.append({"center": list(center), "area_m2": area})
            return centers
        except ImportError:
            return []

    def _find_planar_holes(self, triangles):
        """Scans X-min, X-max, Y-min... faces. If a face has geometry but 'internal' loops, those are vents."""
        all_pts = triangles.reshape(-1, 3)
        min_box = np.min(all_pts, axis=0)
        max_box = np.max(all_pts, axis=0)
        tol = 0.01

        planes = [
            (0, min_box[0]), (0, max_box[0]),
            (1, min_box[1]), (1, max_box[1]),
            (2, min_box[2]), (2, max_box[2])
        ]

        vents = []
        for axis, value in planes:
            mask = np.all(np.abs(triangles[:, :, axis] - value) < tol, axis=1)
            plane_tris = triangles[mask]
            if len(plane_tris) == 0: continue

            edges = defaultdict(int)
            for tri in plane_tris:
                t_2d = [tuple(np.delete(p, axis)) for p in tri]
                t_2d = [tuple(np.round(p, 4)) for p in t_2d]
                pairs = [(t_2d[0], t_2d[1]), (t_2d[1], t_2d[2]), (t_2d[2], t_2d[0])]
                for p1, p2 in pairs:
                    edge = tuple(sorted((p1, p2)))
                    edges[edge] += 1

            boundary_edges = [e for e, count in edges.items() if count == 1]
            if not boundary_edges: continue

            loops = self._order_edges_into_loops(boundary_edges)
            if len(loops) > 1:
                loop_stats = []
                for loop in loops:
                    pts = np.array(loop)
                    mins = np.min(pts, axis=0)
                    maxs = np.max(pts, axis=0)
                    area = np.prod(maxs - mins)
                    center_2d = np.mean(pts, axis=0)

                    center_3d = [0,0,0]
                    center_3d[axis] = value
                    other_axes = [i for i in range(3) if i != axis]
                    center_3d[other_axes[0]] = center_2d[0]
                    center_3d[other_axes[1]] = center_2d[1]

                    loop_stats.append({'area': area, 'center': center_3d})

                loop_stats.sort(key=lambda x: x['area'], reverse=True)
                for vent in loop_stats[1:]:
                    # FIXED: Pack area_m2 into dictionary
                    vents.append({"center": vent['center'], "area_m2": vent['area']})
        return vents

    def _order_edges_into_loops(self, edges):
        adj = defaultdict(list)
        for p1, p2 in edges:
            adj[p1].append(p2)
            adj[p2].append(p1)

        visited_edges = set()
        loops = []

        for start_edge in edges:
            if start_edge in visited_edges: continue
            curr = start_edge[0]
            prev = None
            loop_pts = []

            while True:
                loop_pts.append(curr)
                neighbors = adj[curr]
                next_node = None
                for n in neighbors:
                    edge_key = tuple(sorted((curr, n)))
                    if n != prev and edge_key not in visited_edges:
                        next_node = n
                        visited_edges.add(edge_key)
                        break

                if next_node:
                    prev = curr
                    curr = next_node
                    if curr == loop_pts[0]: break
                else:
                    break

            if len(loop_pts) > 2:
                loops.append(loop_pts)

        return loops

    def _load_stl(self, stl_path):
        triangles = []
        try:
            with open(stl_path, 'r') as f:
                lines = f.readlines()
                if "solid" in lines[0] or "facet" in lines[1]:
                    current = []
                    for line in lines:
                        if 'vertex' in line:
                            parts = line.strip().split()
                            current.append([float(parts[1]), float(parts[2]), float(parts[3])])
                            if len(current) == 3:
                                triangles.append(current)
                                current = []
                    return triangles
        except: pass

        try:
            with open(stl_path, 'rb') as f:
                header = f.read(80)
                count_data = f.read(4)
                count = struct.unpack('<I', count_data)[0]
                for _ in range(count):
                    data = f.read(50)
                    floats = struct.unpack('<12f', data[:48])
                    triangles.append([
                        [floats[3], floats[4], floats[5]],
                        [floats[6], floats[7], floats[8]],
                        [floats[9], floats[10], floats[11]]
                    ])
            return triangles
        except: return []