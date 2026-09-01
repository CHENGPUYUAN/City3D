"""Count open boundary edges in a City3D output mesh (hole detector).

Usage: python3 tools/boundary_check.py <mesh.obj> [more.obj ...]

Reports, per mesh, the number of edges used by exactly one polygon face and
their total length. Point-to-mesh distance metrics barely notice cracks and
holes; a jump in boundary length against a reference output does. The
absolute value is never zero (the extruded model keeps an open rim), so
compare against a known-good output of the same building.
"""
import sys
from collections import Counter


def load_obj(path):
    vertices, faces = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("v "):
                vertices.append(tuple(float(x) for x in line.split()[1:4]))
            elif line.startswith("f "):
                faces.append([int(tok.split("/")[0]) - 1 for tok in line.split()[1:]])
    return vertices, faces


def main() -> None:
    for path in sys.argv[1:]:
        vertices, faces = load_obj(path)
        edge_count = Counter()
        for face in faces:
            for i in range(len(face)):
                edge_count[tuple(sorted((face[i], face[(i + 1) % len(face)])))] += 1
        boundary = [e for e, c in edge_count.items() if c == 1]
        length = sum(
            sum((vertices[a][k] - vertices[b][k]) ** 2 for k in range(3)) ** 0.5
            for a, b in boundary
        )
        print(
            f"{path}: faces={len(faces)} boundary_edges={len(boundary)} "
            f"boundary_len={length:.1f}m"
        )


if __name__ == "__main__":
    main()
