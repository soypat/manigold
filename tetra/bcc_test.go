package tetra

import (
	"testing"

	"gonum.org/v1/gonum/spatial/r3"
)

func TestBCCTetraOrientation(t *testing.T) {
	bcc := MakeBCC(r3.Box{Max: r3.Vec{X: 2, Y: 1, Z: 1}}, 1)
	nodes, tets := bcc.MeshTetraBCC()
	for _, element := range tets {
		tet := [4]r3.Vec{nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]]}
		// tet[1], tet[2] = tet[2], tet[1]
		if !tetraOK(tet) {
			t.Error("X tetrahedron meshed in incorrect order")
		}
	}
	bcc = MakeBCC(r3.Box{Max: r3.Vec{X: 1, Y: 2, Z: 1}}, 1)
	nodes, tets = bcc.MeshTetraBCC()
	for _, element := range tets {
		tet := [4]r3.Vec{nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]]}
		if !tetraOK(tet) {
			t.Error("Y tetrahedron meshed in incorrect order")
		}
	}
	bcc = MakeBCC(r3.Box{Max: r3.Vec{X: 1, Y: 1, Z: 2}}, 1)
	nodes, tets = bcc.MeshTetraBCC()
	for _, element := range tets {
		tet := [4]r3.Vec{nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]]}
		if !tetraOK(tet) {
			t.Error("Z tetrahedron meshed in incorrect order")
		}
	}
}

func tetraOK(t [4]r3.Vec) bool {
	v1 := r3.Sub(t[1], t[0])
	v2 := r3.Sub(t[2], t[0])
	n := r3.Cross(v1, v2)
	v3 := r3.Sub(t[3], t[0])
	return r3.Dot(n, v3) > 0
}
