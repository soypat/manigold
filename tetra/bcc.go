package tetra

import (
	"math"

	"gonum.org/v1/gonum/spatial/r3"
)

// bcc constructs a body centered cubic mesh for isotropic tetrahedron generation.
// Inspired by Tetrahedral Mesh Generation for Deformable Bodies
// Molino, Bridson, Fedkiw.
type bcc struct {
	matrix     bccMatrix
	resolution float64
}

type bccMatrix struct {
	cells []bccCell
	div   [3]int
}

type bccidx int

// BCC node indices. follow same ordering as BoundingBox.Vertices.
const (
	i000 bccidx = iota
	ix00
	ixy0
	i0y0
	i00z
	ix0z
	ixyz
	i0yz
	ictr // BCC central node index.
	nBCC // number of BCC nodes.
)

var unmeshed = [nBCC]int{-1, -1, -1 /**/, -1, -1, -1 /**/, -1, -1, -1}

type bccCell struct {
	bccnod [nBCC]int
	pos    r3.Vec
	// BCC nodes indices on
	parent *bccCell
	xp     *bccCell
	xm     *bccCell
	yp     *bccCell
	ym     *bccCell
	zp     *bccCell
	zm     *bccCell
	m      *bcc
}

func (c *bccCell) nodeAt(idx bccidx) int {
	if c == nil {
		return -1
	}
	if idx >= nBCC {
		panic("bad bcc node index")
	}
	return c.bccnod[idx]
}

// neighborNode returns the global node index at n's local bccidx position
// by checking other neighboring cells which share the node. It returns -1 if node
// does not exist or if it is the central BCC node.
func (c *bccCell) neighborNode(idx bccidx) int {
	var cx, cy, cz int
	switch idx {
	case ictr:
		// central node has no junction.
		return -1
	case i000:
		cx = c.xm.nodeAt(ix00)
		cy = c.ym.nodeAt(i0y0)
		cz = c.zm.nodeAt(i00z)
	case ix00:
		cx = c.xp.nodeAt(i000)
		cy = c.ym.nodeAt(ixy0)
		cz = c.zm.nodeAt(ix0z)
	case ixy0:
		cx = c.xp.nodeAt(i0y0)
		cy = c.yp.nodeAt(ix00)
		cz = c.zm.nodeAt(ixyz)
	case i0y0:
		cx = c.xm.nodeAt(ixy0)
		cy = c.yp.nodeAt(i000)
		cz = c.zm.nodeAt(i0yz)
	case i00z:
		cx = c.xm.nodeAt(ix0z)
		cy = c.ym.nodeAt(i0yz)
		cz = c.zp.nodeAt(i000)
	case ix0z:
		cx = c.xp.nodeAt(i00z)
		cy = c.ym.nodeAt(ixyz)
		cz = c.zp.nodeAt(ix00)
	case ixyz:
		cx = c.xp.nodeAt(i0yz)
		cy = c.yp.nodeAt(ix0z)
		cz = c.zp.nodeAt(ixy0)
	case i0yz:
		cx = c.xm.nodeAt(ixyz)
		cy = c.yp.nodeAt(i00z)
		cz = c.zp.nodeAt(i0y0)
	}
	bad := cx >= 0 && cy >= 0 && cx != cy ||
		cx >= 0 && cz >= 0 && cx != cz ||
		cz >= 0 && cy >= 0 && cz != cy
	if bad {
		panic("bad mesh operation detected")
	}
	return max(cx, max(cy, cz))
}

func MakeBCC(b r3.Box, resolution float64) *bcc {
	sz := boxSize(b)
	bccDiv := [3]int{
		int(math.Round(sz.X / resolution)),
		int(math.Round(sz.Y / resolution)),
		int(math.Round(sz.Z / resolution)),
	}
	if bccDiv[0] < 1 || bccDiv[1] < 1 || bccDiv[2] < 1 {
		panic("resolution too low")
	}
	// Number of BCC cells
	Ncells := bccDiv[0] * bccDiv[1] * bccDiv[2]
	mesh := &bcc{
		resolution: resolution,
	}
	matrix := bccMatrix{cells: make([]bccCell, Ncells), div: bccDiv}
	for i := 0; i < bccDiv[0]; i++ {
		x := (float64(i)+0.5)*resolution + b.Min.X
		for j := 0; j < bccDiv[1]; j++ {
			y := (float64(j)+0.5)*resolution + b.Min.Y
			for k := 0; k < bccDiv[2]; k++ {
				z := (float64(k)+0.5)*resolution + b.Min.Z
				matrix.setLevel0(i, j, k, bccCell{pos: r3.Vec{X: x, Y: y, Z: z}, m: mesh, bccnod: unmeshed})
			}
		}
	}
	mesh.matrix = matrix
	return mesh
}

func (t *bcc) MeshTetraBCC() (nodes []r3.Vec, tetras [][4]int) {
	n := 0
	tetras = make([][4]int, 0, len(t.matrix.cells))
	t.matrix.foreach(func(_, _, _ int, cell *bccCell) {
		bb := cell.box()
		vert := boxVertices(bb)
		nctr := n
		cell.bccnod[ictr] = nctr
		n++
		nodes = append(nodes, boxCenter(bb))
		for in := i000; in < ictr; in++ {
			v := cell.neighborNode(in)
			if v == -1 {
				// If node has not yet been created we create it here.
				cell.bccnod[in] = n
				n++
				nodes = append(nodes, vert[in])
			} else {
				cell.bccnod[in] = v
			}
		}
		tetras = append(tetras, cell.tetras()...)
	})
	return nodes, tetras
}

// exists returns true if bccCell is initialized and exists in a mesh.
// Returns false if called on nil bccCell.
func (t *bccCell) exists() bool {
	return t != nil && t.m != nil
}

func (t *bccCell) box() r3.Box {
	res := t.m.resolution
	return centeredBox(t.pos, r3.Vec{X: res, Y: res, Z: res})
}

func (m *bccMatrix) setLevel0(i, j, k int, c bccCell) {
	if i < 0 || j < 0 || k < 0 || i >= m.div[0] || j >= m.div[1] || k >= m.div[2] {
		panic("oob tmatrix access")
	}
	ca := m.at(i, j, k)
	*ca = c
	// Update x neighbors
	ca.xm = m.at(i-1, j, k)
	if ca.xm.exists() {
		ca.xm.xp = ca
	}
	ca.xp = m.at(i+1, j, k)
	if ca.xp.exists() {
		ca.xp.xm = ca
	}
	// Update y neighbors.
	ca.ym = m.at(i, j-1, k)
	if ca.ym.exists() {
		ca.ym.yp = ca
	}
	ca.yp = m.at(i, j+1, k)
	if ca.yp.exists() {
		ca.yp.ym = ca
	}
	// Update z neighbors
	ca.zm = m.at(i, j, k-1)
	if ca.zm.exists() {
		ca.zm.zp = ca
	}
	ca.zp = m.at(i, j, k+1)
	if ca.zp.exists() {
		ca.zp.zm = ca
	}
}

func (m *bccMatrix) at(i, j, k int) *bccCell {
	if i < 0 || j < 0 || k < 0 || i >= m.div[0] || j >= m.div[1] || k >= m.div[2] {
		return nil
	}
	return &m.cells[i*m.div[1]*m.div[2]+j*m.div[2]+k]
}

func (m *bccMatrix) foreach(f func(i, j, k int, nod *bccCell)) {
	for i := 0; i < m.div[0]; i++ {
		ii := i * m.div[1] * m.div[2]
		for j := 0; j < m.div[1]; j++ {
			jj := j * m.div[2]
			for k := 0; k < m.div[2]; k++ {
				f(i, j, k, &m.cells[ii+jj+k])
			}
		}
	}
}

func max(a, b int) int {
	if a >= b {
		return a
	}
	return b
}

// naive implementation of BCC tetra mesher. Meshing is not isotropic.
func (c bccCell) naiveTetras() [][4]int {
	nctr := c.bccnod[ictr]
	return [][4]int{
		// YZ plane facing tetrahedrons.
		{c.bccnod[i000], c.bccnod[i0yz], c.bccnod[i0y0], nctr},
		{c.bccnod[i000], c.bccnod[i00z], c.bccnod[i0yz], nctr},
		{c.bccnod[ix00], c.bccnod[ixy0], c.bccnod[ixyz], nctr},
		{c.bccnod[ix00], c.bccnod[ixyz], c.bccnod[ix0z], nctr},
		// XZ
		{c.bccnod[i000], c.bccnod[ix0z], c.bccnod[i00z], nctr},
		{c.bccnod[i000], c.bccnod[ix00], c.bccnod[ix0z], nctr},
		{c.bccnod[i0y0], c.bccnod[i0yz], c.bccnod[ixyz], nctr},
		{c.bccnod[i0y0], c.bccnod[ixyz], c.bccnod[ixy0], nctr},
		// XY
		{c.bccnod[i000], c.bccnod[ixy0], c.bccnod[ix00], nctr},
		{c.bccnod[i000], c.bccnod[i0y0], c.bccnod[ixy0], nctr},
		{c.bccnod[i00z], c.bccnod[ix0z], c.bccnod[ixyz], nctr},
		{c.bccnod[i00z], c.bccnod[ixyz], c.bccnod[i0yz], nctr},
	}
}

// tetras is the BCC lattice meshing method. Results in isotropic mesh.
func (c bccCell) tetras() (tetras [][4]int) {
	// We mesh tetrahedrons on minor sides.
	nctr := c.bccnod[ictr]
	// Start with nodes in z direction since matrix is indexed with z as major
	// dimension so maybe zm is on the cache.
	if c.zm.exists() && c.zm.bccnod[ictr] >= 0 {
		zctr := c.zm.bccnod[ictr]
		tetras = append(tetras,
			[4]int{nctr, c.bccnod[i000], c.bccnod[ix00], zctr},
			[4]int{nctr, c.bccnod[ix00], c.bccnod[ixy0], zctr},
			[4]int{nctr, c.bccnod[ixy0], c.bccnod[i0y0], zctr},
			[4]int{nctr, c.bccnod[i0y0], c.bccnod[i000], zctr},
		)
	}
	if c.ym.exists() && c.ym.bccnod[ictr] >= 0 {
		yctr := c.ym.bccnod[ictr]
		tetras = append(tetras,
			[4]int{nctr, c.bccnod[ix00], c.bccnod[i000], yctr},
			[4]int{nctr, c.bccnod[ix0z], c.bccnod[ix00], yctr},
			[4]int{nctr, c.bccnod[i00z], c.bccnod[ix0z], yctr},
			[4]int{nctr, c.bccnod[i000], c.bccnod[i00z], yctr},
		)
	}
	if c.xm.exists() && c.xm.bccnod[ictr] >= 0 {
		xctr := c.xm.bccnod[ictr]
		tetras = append(tetras,
			[4]int{nctr, c.bccnod[i000], c.bccnod[i0y0], xctr},
			[4]int{nctr, c.bccnod[i00z], c.bccnod[i000], xctr},
			[4]int{nctr, c.bccnod[i0yz], c.bccnod[i00z], xctr},
			[4]int{nctr, c.bccnod[i0y0], c.bccnod[i0yz], xctr},
		)
	}
	return tetras
}
