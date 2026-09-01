"""Point-to-mesh distance between an input cloud and a City3D output mesh.

Usage: uv run --python 3.12 --with open3d --with mapbox-earcut tools/eval_quality.py <cloud.ply> <mesh.obj> <footprint.obj>

Prints the distance from every input point inside the footprint polygon(s) to
the closest mesh triangle (median / mean / p90 / p95 in meters). All files
are expected in the same local coordinate frame (the CLI writes meshes
without the offset shift). The OBJ files are parsed here (open3d's reader
skips the polygon faces City3D writes); polygon faces are triangulated with
earcut — City3D outputs concave n-gons, which a fan triangulation would
corrupt.
"""
import sys

import numpy as np
import open3d as o3d
from mapbox_earcut import triangulate_float32


def load_obj(path: str) -> tuple[list[tuple[float, float, float]], list[list[int]]]:
    vertices: list[tuple[float, float, float]] = []
    faces: list[list[int]] = []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                _, x, y, z = line.split()[:4]
                vertices.append((float(x), float(y), float(z)))
            elif line.startswith("f "):
                faces.append([int(tok.split("/")[0]) - 1 for tok in line.split()[1:]])
    return vertices, faces


def polygon_triangles(face: list[int], verts: np.ndarray) -> list[tuple[int, int, int]]:
    if len(face) == 3:
        return [tuple(face)]
    ring = verts[face]
    # project the ring onto its dominant plane for a 2D earcut
    normal = np.cross(ring[1] - ring[0], ring[-1] - ring[0])
    axis = np.argsort(np.abs(normal))
    u, v = axis[0], axis[1]  # the two axes most perpendicular to the normal
    flat = np.ascontiguousarray(np.array([ring[:, u], ring[:, v]], dtype=np.float32).T)
    tri_idx = triangulate_float32(flat, np.array([len(face)], dtype=np.uint32))
    return [
        (face[tri_idx[i]], face[tri_idx[i + 1]], face[tri_idx[i + 2]])
        for i in range(0, len(tri_idx), 3)
    ]


def points_in_polygon(pts_xy: np.ndarray, poly_xy: np.ndarray) -> np.ndarray:
    """Even-odd rule, vectorized over points."""
    x, y = pts_xy[:, 0], pts_xy[:, 1]
    inside = np.zeros(len(pts_xy), dtype=bool)
    for i in range(len(poly_xy)):
        x1, y1 = poly_xy[i]
        x2, y2 = poly_xy[(i + 1) % len(poly_xy)]
        crosses = (y1 > y) != (y2 > y)
        with np.errstate(divide="ignore", invalid="ignore"):
            xint = np.where(crosses, (x2 - x1) * (y - y1) / (y2 - y1) + x1, np.inf)
        inside ^= crosses & (x < xint)
    return inside


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    cloud_file, mesh_file, footprint_file = sys.argv[1:4]

    pcd = o3d.io.read_point_cloud(cloud_file)
    verts, mesh_faces = load_obj(mesh_file)
    fp_verts, fp_faces = load_obj(footprint_file)
    verts_arr = np.asarray(verts)
    tris = [t for f in mesh_faces for t in polygon_triangles(f, verts_arr)]
    if len(pcd.points) == 0 or len(tris) == 0:
        sys.exit(f"empty input: points={len(pcd.points)} triangles={len(tris)}")

    pts = np.asarray(pcd.points)
    keep = np.zeros(len(pts), dtype=bool)
    for f in fp_faces:
        poly = np.asarray([[fp_verts[i][0], fp_verts[i][1]] for i in f])
        keep |= points_in_polygon(pts[:, :2], poly)
    pts = pts[keep]

    mesh = o3d.t.geometry.TriangleMesh(
        o3d.core.Tensor(np.asarray(verts, dtype=np.float32)),
        o3d.core.Tensor(np.asarray(tris, dtype=np.int32)),
    )
    scene = o3d.t.geometry.RaycastingScene()
    _ = scene.add_triangles(mesh)
    d = scene.compute_distance(o3d.core.Tensor(pts.astype(np.float32))).numpy()

    print(
        f"points={len(d)}/{keep.sum()} triangles={len(tris)} "
        f"median={np.median(d):.4f} mean={d.mean():.4f} "
        f"p90={np.percentile(d, 90):.4f} p95={np.percentile(d, 95):.4f} max={d.max():.4f}"
    )


if __name__ == "__main__":
    main()
