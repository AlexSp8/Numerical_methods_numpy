
from numpy.typing import NDArray
import numpy as np

class MeshDiscretization():

    def __init__(self, x0: NDArray[np.float64], xf: NDArray[np.float64]):

        self.x0 = np.array(x0, dtype=float)
        self.xf = np.array(xf, dtype=float)
        self.ndim = x0.shape[0]
        self.nodes_total = 0
        self.x_mesh = None
        self.connectivity = None

    def set_rectangular_mesh_coordinates(self, nodes_dim: NDArray[np.int32],
        p: NDArray[np.float64] = None, major_order: str = 'row'):
        """Sets the mesh coordinates for a rectangular mesh
        given the stretching factor for refinement in each direction
        and the order of node enumeration.

        Args:
            self: the MeshDiscretization object that holds (x0, xf, ndim)
            nodes_dim: the number of nodes in each direction
            p (optional): the power law refinement factor in each direction.
                For p > 1: refinement towards x0.
                For 0 < p < 1: refinement towards xf
            major_order (optional): the order of enumeration
        Returns:
            The mesh coordinates array in self.x_mesh."""

        ndim = self.ndim

        if p is None:
            p = np.ones(ndim, dtype=float)

        nodes_total = 1
        for d in range(ndim):
            nodes_total *= nodes_dim[d]

        self.nodes_total = int(nodes_total)

        self.x_mesh = np.zeros((nodes_total,ndim))

        x_1d = []
        for d in range(ndim):
            xd = np.array([self.x0[d], self.xf[d]])
            nd = nodes_dim[d]
            x_1d.append(self.create_mesh_1d(nd, xd, p[d]))

        if major_order == 'row':
            # row-major order: the first dimension (x) cycles first
            dim_start, dim_end, dim_step = 0, ndim, +1
        else:
            # column-major order: the last dimension (z) cycles first
            dim_start, dim_end, dim_step = ndim-1, -1, -1

        for inod in range(nodes_total):
            idx = inod
            for d in range(dim_start, dim_end, dim_step):
                nd = nodes_dim[d]
                idx_dim = idx%nd
                self.x_mesh[inod,d] = x_1d[d][idx_dim]
                idx = idx//nd

    @staticmethod
    def create_mesh_1d(nodes: int, xd: NDArray[np.float64],
        p: int = 1) -> NDArray[np.float64]:
        """Returns the mesh coordinates given the number of nodes,
        the domain bounds and the stretching factor for refinement.

        Args:
            nodes: the number of nodes
            xd: the domain bounds
            p (optional): the power law refinement factor.
                For p > 1: refinement towards x0.
                For 0 < p < 1: refinement towards xf
        Returns:
            The array of 1D mesh coordinates."""

        x0, xf = xd[0], xd[1]
        L = xf - x0

        x = np.zeros(nodes)
        for i in range(nodes):
            xi = i/(nodes-1)
            gi = xi**p
            # gi = 0.5*(1.0-np.cos(p*np.pi*xi))
            # gi = np.tanh(p*(2*xi-1))/(2*np.tanh(p)) + 0.5
            x[i] = x0 + L*gi

        return x
        # xi = np.linspace(0, 1, nodes)
        # gi = xi**p
        # return x0 + L * gi

    def check_dt(self, D0: float, v0: float, dt0: float) -> float:
        """Checks if the time discretization is stable.

        Args:
            x: the spatial domain
            D0: the diffusion number
            v0: the convection velocity
            dt0: the initial time step to check
        Returns:
            The initial time step if it is stable. Otherwise the maximum stable dt."""

        dt = dt0

        D = np.maximum(np.abs(D0), 1e-12)
        v = np.maximum(np.abs(v0), 1e-12)

        # dx = np.diff(self.x_mesh)
        dx = self.x_mesh[1:] - self.x_mesh[0:-1]
        dx_min = np.min(dx)
        dx_max = np.max(dx)

        dt_max_diffusion = (dx_min**2)/(2*D)
        dt_max_convection = (2*D)/(v**2)
        dt_max = min(dt_max_diffusion, dt_max_convection)

        cfl_diffusion = dt/(dx_min**2)
        courant = (v*dt)/dx_min

        print("="*40)
        print("    CONVECTION-DIFFUSION DIAGNOSTICS")
        print("="*40)
        print(f"dx_min:              {dx_min:.4f}")
        print(f"dx_max:              {dx_max:.4f}\n")
        print(f"Courant (C):         {courant:.4f}\n")
        print(f"Max dt (Diffusion):  {dt_max_diffusion:.6f}")
        print(f"Max dt (Convection): {dt_max_convection:.6f}")
        print(f"Max dt Allowed:      {dt_max:.6f}")
        print(f"Time Step (dt):      {dt:.6f}")
        print("-"*40)

        if dt > dt_max:
            print(f"dt = {dt:.6f} exceeds the maximum stable dt_max = {dt_max:.6f}.")
            print(f"Setting dt = {dt_max:.6f}.")
            dt = dt_max
            courant = (v*dt)/dx_min
            print(f"Courant (C):         {courant:.4f}\n")

        print("="*40)

        return dt

    def check_dx(self, D0: float, v0: float) -> None:
        """Checks if the spatial discretization is stable.

        Args:
            x: the spatial domain
            D0: the diffusion number
            v0: the convection velocity
        Returns:
            The initial time step if it is stable. Otherwise the maximum stable dt."""

        D = np.maximum(np.abs(D0), 1e-12)
        v = np.maximum(np.abs(v0), 1e-12)

        # dx = np.diff(self.x_mesh)
        dx = self.x_mesh[1:] - self.x_mesh[0:-1]
        dx_max = np.max(dx)

        pe = (v*dx_max)/D

        dx_stable = D*2.0/v
        if pe > 2.0:
            print(f"dx = {dx_max:.4f} > dx_max = {dx_stable:.4f} (Pe = {pe:.1f} > 2.0)")
            raise ValueError('Too large dx for this system')


class MeshDiscretizationFD(MeshDiscretization):

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)
        self.bulk_nodes = []
        self.boundary_nodes = []

    def set_rectangular_mesh_fd(self, nodes_dim: NDArray[np.int32] = None,
        p: NDArray[np.float64] = None, major_order: str = 'row'):

        self.set_rectangular_mesh_coordinates(nodes_dim, p, major_order)
        self.set_rectangular_mesh_connectivity_fd(nodes_dim, major_order)
        self.set_rectangular_mesh_node_types_fd()

    def set_rectangular_mesh_connectivity_fd(self, nodes_dim: NDArray[np.int32] = None,
        major_order: str = 'row'):
        """Sets the connectivity for a FD rectangular mesh given
        the number of nodes in each direction and order of node enumeration.

        Args:
            self: the MeshDiscretization object that holds (ndim, nodes_total)
            nodes_dim: the number of nodes in each direction
            major_order (optional): the order of enumeration
        Returns:
            The FD connectivity array in self.connectivity."""

        ndim = self.ndim
        nodes_total = self.nodes_total

        self.connectivity = np.full((nodes_total, ndim, 2), -1, dtype=int)

        node_step = [1]*ndim
        if major_order == 'row':
            # row-major order: the first dimension (x) cycles first
            for d in range(1, ndim):
                node_step[d] = node_step[d-1]*nodes_dim[d-1]
            dim_start, dim_end, dim_step = 0, ndim, +1
        else:
            # column-major order: the last dimension (z) cycles first
            for d in range(ndim - 2, -1, -1):
                node_step[d] = node_step[d+1]*nodes_dim[d+1]
            dim_start, dim_end, dim_step = ndim-1, -1, -1

        for inod in range(nodes_total):

            idx_dim = [0]*ndim
            idx = inod
            for d in range(dim_start, dim_end, dim_step):
                nd = nodes_dim[d]
                idx_dim[d] = idx%nd
                idx = idx//nd

            for d in range(ndim):
                # Backward Neighbor check
                if idx_dim[d] > 0:
                    self.connectivity[inod, d, 0] = inod - node_step[d]
                # Forward Neighbor check
                if idx_dim[d] < nodes_dim[d] - 1:
                    self.connectivity[inod, d, 1] = inod + node_step[d]

    def set_rectangular_mesh_node_types_fd(self):
        """Classifies the nodes for a FD rectangular mesh (bulk, boundary, multi-boundary).

        Args:
            self: the MeshDiscretization object that holds (nodes_total, ndim, connectivity)
        Returns:
            The lists of node types in self for each type of node
            (bulk_nodes, boundary_nodes)."""

        nodes_total = self.nodes_total
        ndim = self.ndim
        connectivity = self.connectivity

        self.bulk_nodes = []
        self.boundary_nodes = []

        for gnod in range(nodes_total):

            boundaries = []
            for d in range(ndim):
                for j in range(connectivity.shape[2]):
                    if connectivity[gnod, d, j] == -1:
                        bnd = 2*d+j
                        boundaries.append(bnd)
                        break

            bnd_count = len(boundaries)
            if bnd_count == 0:
                self.bulk_nodes.append(gnod)
            else:
                self.boundary_nodes.append((gnod, boundaries))


class MeshDiscretizationFE(MeshDiscretization):

    def __init__(self, nbf: int, fe_type: str = 'line', *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.nbf = nbf
        self.fe_type = fe_type

        self._set_fe_order()

        self._set_fe_edges()
        self._set_edge_nodes()

        self._set_fe_faces()
        self._set_face_nodes()

        self.nel = 0
        self.bulk_nodes = []
        self.boundary_nodes = []

    def _set_fe_order(self):
        """Sets the order of finite elements.

        Args:
            self: the MeshDiscretization object that holds (nbf, fe_type)
        Returns:
            The FE order in self.fe_order."""

        fe_type = self.fe_type
        nbf = self.nbf

        fe_order = None

        if fe_type == 'line':
            if nbf == 2:
                fe_order = 'linear'
            elif nbf == 3:
                fe_order = 'quadratic'
            else:
                raise NotImplementedError(f'Unsupported nbf = {nbf} for line in _set_fe_order!')

        elif fe_type == 'quadrangle':
            if nbf == 4:
                fe_order = 'linear'
            elif nbf == 9:
                fe_order = 'quadratic'
            elif nbf == 8:
                fe_order = 'serendipity'
            else:
                raise NotImplementedError(f'Unsupported nbf = {nbf} for quadrangle in _set_fe_order!')

        elif fe_type == 'triangle':
            if nbf == 3:
                fe_order = 'linear'
            elif nbf == 6:
                fe_order = 'quadratic'
            else:
                raise NotImplementedError(f'Unsupported nbf = {nbf} for triangle in _set_fe_order!')

        elif fe_type == 'hexahedron':
            if nbf == 8:
                fe_order = 'linear'
            elif nbf == 27:
                fe_order = 'quadratic'
            elif nbf == 20:
                fe_order = 'serendipity'
            else:
                raise NotImplementedError(f'Unsupported nbf = {nbf} for hexahedron in _set_fe_order!')

        elif fe_type == 'tetrahedron':
            if nbf == 4:
                fe_order = 'linear'
            elif nbf == 10:
                fe_order = 'quadratic'
            else:
                raise NotImplementedError(f'Unsupported nbf = {nbf} for tetrahedron in _set_fe_order!')
        else:
            raise NotImplementedError(f'Unsupported fe_type = {fe_type} in _set_fe_order!')

        self.fe_order = fe_order

    def _set_fe_edges(self):
        """Sets the number of edges of the element.

        Args:
            self: the MeshDiscretization object that holds (fe_type)
        Returns:
            The number of edges in self.ned."""

        fe_type = self.fe_type

        ned = 0
        if fe_type == 'line':
            ned = 0
        elif fe_type == 'quadrangle':
            ned = 4
        elif fe_type == 'triangle':
            ned = 3
        elif fe_type == 'hexahedron':
            ned = 12
        elif fe_type == 'tetrahedron':
            ned = 6
        else:
            raise NotImplementedError(f'Unsupported fe_type = {fe_type} in _set_fe_edges!')

        self.ned = ned

    def _set_fe_faces(self):
        """Sets the number of faces of the element.

        Args:
            self: the MeshDiscretization object that holds (fe_type)
        Returns:
            The number of faces in self.nfc."""

        fe_type = self.fe_type

        nfc = 0
        if fe_type == 'line':
            nfc = 0
        elif fe_type == 'quadrangle':
            nfc = 0
        elif fe_type == 'triangle':
            nfc = 0
        elif fe_type == 'hexahedron':
            nfc = 6
        elif fe_type == 'tetrahedron':
            nfc = 4
        else:
            raise NotImplementedError(f'Unsupported fe_type = {fe_type} in _set_fe_faces!')

        self.nfc = nfc

    def _set_edge_nodes(self):
        """Sets the node numbering on the element edges.

        Args:
            self: the MeshDiscretization object that holds (fe_type, fe_order, ned)
        Returns:
            The node numbering on each edge of the element."""

        fe_type = self.fe_type
        fe_order = self.fe_order

        if fe_order == 'linear':
            nbf_1d = 2
        elif fe_order == 'quadratic':
            nbf_1d = 3

        edges = np.zeros((self.ned, nbf_1d))

        if fe_type == 'line':
            pass

        elif fe_type == 'quadrangle':
            if fe_order == 'linear':
                edges[0,:] = [1,2] # e = -1
                edges[1,:] = [2,3] # c = +1
                edges[2,:] = [3,4] # e = +1
                edges[3,:] = [4,1] # c = -1
            elif fe_order == 'quadratic' or fe_order == 'serendipity':
                edges[0,:] = [1,2,5] # e = -1
                edges[1,:] = [2,3,6] # c = +1
                edges[2,:] = [3,4,7] # e = +1
                edges[3,:] = [4,1,8] # c = -1

        elif fe_type == 'triangle':
            if fe_order == 'linear':
                edges[0,:] = [1,2] # e = 0
                edges[1,:] = [1,3] # c = 0
                edges[2,:] = [2,3] # e = 1-c
            elif fe_order == 'quadratic' or fe_order == 'serendipity':
                edges[0,:] = [1,2,4] # e = 0
                edges[1,:] = [1,3,6] # c = 0
                edges[2,:] = [2,3,5] # e = 1-c

        elif fe_type == 'hexahedron':
            if fe_order == 'linear':
                edges[0,:] = [1,2] # e = -1, s = -1
                edges[1,:] = [2,3] # c = +1, s = -1
                edges[2,:] = [3,4] # e = +1, s = -1
                edges[3,:] = [4,1] # c = -1, s = -1
                edges[4,:] = [5,6] # e = -1, s = +1
                edges[5,:] = [6,7] # c = +1, s = +1
                edges[6,:] = [7,8] # e = +1, s = +1
                edges[7,:] = [8,5] # c = -1, s = +1
                edges[8,:] = [1,5] # c = -1, e = -1
                edges[9,:] = [2,6] # c = +1, e = -1
                edges[10,:] = [3,7] # c = +1, e = +1
                edges[11,:] = [4,8] # c = -1, e = +1
            elif fe_order == 'quadratic' or fe_order == 'serendipity':
                edges[0,:] = [1,2,9] # e = -1, s = -1
                edges[1,:] = [2,3,10] # c = +1, s = -1
                edges[2,:] = [3,4,11] # e = +1, s = -1
                edges[3,:] = [4,1,12] # c = -1, s = -1
                edges[4,:] = [5,6,13] # e = -1, s = +1
                edges[5,:] = [6,7,14] # c = +1, s = +1
                edges[6,:] = [7,8,15] # e = +1, s = +1
                edges[7,:] = [8,5,16] # c = -1, s = +1
                edges[8,:] = [1,5,17] # c = -1, e = -1
                edges[9,:] = [2,6,18] # c = +1, e = -1
                edges[10,:] = [3,7,19] # c = +1, e = +1
                edges[11,:] = [4,8,20] # c = -1, e = +1

        elif fe_type == 'tetrahedron':
            if fe_order == 'linear':
                edges[0,:] = [1,2] # e = 0, s = 0
                edges[1,:] = [1,3] # c = 0, s = 0
                edges[2,:] = [1,4] # c = 0, e = 0
                edges[3,:] = [2,3] # e = 1-c, s = 0
                edges[4,:] = [3,4] # c = 0, s = 1-e
                edges[5,:] = [2,4] # e = 0, s = 1-c
            elif fe_order == 'quadratic' or fe_order == 'serendipity':
                edges[0,:] = [1,2,5] # e = 0, s = 0
                edges[1,:] = [1,3,7] # c = 0, s = 0
                edges[2,:] = [1,4,8] # c = 0, e = 0
                edges[3,:] = [2,3,6] # e = 1-c, s = 0
                edges[4,:] = [3,4,10] # c = 0, s = 1-e
                edges[5,:] = [2,4,9] # e = 0, s = 1-c
        else:
            raise NotImplementedError(f'Unsupported fe_type = {fe_type} in _set_edge_nodes!')

        self.edges = edges

    def _set_face_nodes(self):
        """Sets the node numbering on the element faces.

        Args:
            self: the MeshDiscretization object that holds (fe_type, fe_order, ned)
        Returns:
            The node numbering on each face of the element."""

        fe_type = self.fe_type
        fe_order = self.fe_order

        nfc = self.nfc

        if fe_type == 'line':
            faces = np.zeros((0, 0))

        elif fe_type == 'quadrangle':
            faces = np.zeros((0, 0))

        elif fe_type == 'triangle':
            faces = np.zeros((0, 0))

        elif fe_type == 'hexahedron':
            if fe_order == 'linear':
                nbf_2d = 4
                faces = np.zeros((nfc, nbf_2d))
                faces[0,:] = [1,2,3,4] # s = -1
                faces[1,:] = [5,6,7,8] # s = +1
                faces[2,:] = [1,2,6,5] # e = -1
                faces[3,:] = [4,3,7,8] # e = +1
                faces[4,:] = [1,4,8,5] # c = -1
                faces[5,:] = [2,3,7,6] # c = +1
            elif fe_order == 'quadratic':
                nbf_2d = 9
                faces = np.zeros((nfc, nbf_2d))
                faces[0,:] = [1,2,3,4,9,10,11,12,25] # s = -1
                # faces[0,:] = [1,2,3,4,9,10,11,12,21] # s = -1, Salome
                faces[1,:] = [5,6,7,8,13,14,15,16,26] # s = +1
                faces[2,:] = [1,2,6,5,9,18,13,17,21] # e = -1
                # faces[2,:] = [1,2,6,5,9,18,13,17,22] # e = -1, Salome
                faces[3,:] = [4,3,7,8,11,19,15,20,23] # e = +1
                # faces[3,:] = [4,3,7,8,11,19,15,20,24] # e = +1, Salome
                faces[4,:] = [1,4,8,5,12,20,16,17,24] # c = -1
                # faces[4,:] = [1,4,8,5,12,20,16,17,25] # c = -1, Salome
                faces[5,:] = [2,3,7,6,10,19,14,18,22] # c = +1
                # faces[5,:] = [2,3,7,6,10,19,14,18,23] # c = +1, Salome
            elif fe_order == 'serendipity':
                nbf_2d = 8
                faces = np.zeros((nfc, nbf_2d))
                faces[0,:] = [1,2,3,4,9,10,11,12] # s = -1
                faces[1,:] = [5,6,7,8,13,14,15,16] # s = +1
                faces[2,:] = [1,2,6,5,9,18,13,17] # e = -1
                faces[3,:] = [4,3,7,8,11,19,15,20] # e = +1
                faces[4,:] = [1,4,8,5,12,20,16,17] # c = -1
                faces[5,:] = [2,3,7,6,10,19,14,18] # c = +1

        elif fe_type == 'tetrahedron':
            if fe_order == 'linear':
                nbf_2d = 3
                faces = np.zeros((nfc, nbf_2d))
                faces[0,:] = [1,3,4]
                faces[1,:] = [1,4,2]
                faces[2,:] = [1,2,3]
                faces[3,:] = [2,3,4]
            elif fe_order == 'quadratic' or fe_order == 'serendipity':
                nbf_2d = 6
                faces = np.zeros((nfc, nbf_2d))
                faces[0,:] = [1,3,4,7,10,8]
                faces[1,:] = [1,4,2,5,9,8]
                faces[2,:] = [1,2,3,5,6,7]
                faces[3,:] = [2,3,4,6,9,10]
        else:
            raise NotImplementedError(f'Unsupported fe_type = {fe_type} in _set_face_nodes!')

        self.faces = faces

    def set_rectangular_mesh_fe(self, nel_dim: NDArray[np.int32],
        p: NDArray[np.float64] = None, major_order: str = 'row'):
        """Sets up a rectangular structured FE mesh.

        Args:
            self: the MeshDiscretization object that holds (nbf, fe_type)
            nel_dim: the array of the number of elements in each dimension
            p (optional): the power law refinement factor in each direction.
                For p > 1: refinement towards x0.
                For 0 < p < 1: refinement towards xf
            major_order (optional): the order of enumeration
        Returns:
            The structure FE mesh data in self (nel, x_mesh, connectivity,
            boundary_nodes)."""

        fe_type = self.fe_type
        fe_order = self.fe_order

        if fe_type in ['line', 'quadrangle', 'hexahedron']:
            nel_per_cell = 1
        elif fe_type == 'triangle':
            nel_per_cell = 2
        elif fe_type == 'tetrahedron':
            nel_per_cell = 6
        else:
            raise NotImplementedError(
                f"Unsupported fe_type: {fe_type} in set_rectangular_mesh_fe!")

        self.nel = 1
        for d in range(self.ndim):
            self.nel *= nel_dim[d]

        self.nel *= nel_per_cell

        if fe_order == 'linear':
            nbf_1d = 2
        elif fe_order in ['quadratic', 'serendipity']:
            nbf_1d = 3

        nodes_dim = (nbf_1d-1)*nel_dim + 1
        self.set_rectangular_mesh_coordinates(nodes_dim, p, major_order)

        if fe_order == 'serendipity':
            self.set_rectangular_mesh_serendipity(nel_dim, nel_per_cell, major_order)
        else:
            self.set_rectangular_mesh_connectivity_fe(nel_dim, nel_per_cell, major_order)

        self.set_rectangular_mesh_node_types_fe()

    def set_rectangular_mesh_connectivity_fe(self, nel_dim: NDArray[np.int32],
        nel_per_cell: int, major_order: str = 'row'):
        """Sets the connectivity for a FE rectangular mesh
        given the order of node enumeration.

        Args:
            self: the MeshDiscretization object that holds (nbf, nel)
            nel_dim: the array of the number of elements in each dimension
            nel_per_cell: the number of elements in which a cell is divided
            major_order (optional): the order of enumeration
        Returns:
            The connectivity array in self.connectivity."""

        fe_order = self.fe_order
        ndim = self.ndim
        nbf = self.nbf
        nel = self.nel

        if fe_order == 'linear':
            nbf_1d = 2
            node_step = 1
        elif fe_order in ['quadratic', 'serendipity']:
            nbf_1d = 3
            node_step = 2
        else:
            raise NotImplementedError(
                f"Unsupported fe_order = {fe_order} in set_rectangular_mesh_connectivity_fe!")

        self.connectivity = np.zeros((nel, nbf), dtype=int)

        nodes_dim = (nbf_1d-1)*nel_dim + 1

        local_nodes = self.set_local_node_indexing(nel_per_cell)

        if major_order == 'row':
            # row-major order: the first dimension (x) cycles first
            dim_start, dim_end, dim_step = 0, ndim, +1
        else:
            # column-major order: the last dimension (z) cycles first
            dim_start, dim_end, dim_step = ndim-1, -1, -1

        for iel in range(nel):

            gnodes = np.zeros(nbf, dtype=int)
            nodes_product = 1

            remainder = iel
            icell = remainder%nel_per_cell
            remainder //= nel_per_cell

            for d in range(dim_start, dim_end, dim_step):

                nel_d = nel_dim[d]
                iel_idx_dim = remainder%nel_d
                remainder //= nel_d

                nodes_idx = local_nodes[icell,d,:] + iel_idx_dim*node_step
                gnodes += nodes_idx*nodes_product
                nodes_product *= nodes_dim[d]

            self.connectivity[iel,:] = gnodes

    def set_local_node_indexing(self, nel_per_cell: int):
        """Sets the local node indexing for the given element type.

        Args:
            self: the MeshDiscretization object that holds (nbf, nel)
            nel_per_cell: the number of elements in which a cell is divided
        Returns:
            The local node indexing inside the element."""

        ndim = self.ndim
        nbf = self.nbf
        fe_order = self.fe_order
        fe_type = self.fe_type

        local_nodes = np.zeros((nel_per_cell,ndim,nbf), dtype=int)

        if ndim == 1:

            if fe_order == 'linear':
                local_nodes[0,0,:] = [0,1]
            elif fe_order == 'quadratic':
                local_nodes[0,0,:] = [0,1,2]
            else:
                raise NotImplementedError(
                    f"Unsupported fe_order: {fe_order} for ndim = 1 in set_local_node_indexing!")

        elif ndim == 2:

            if fe_type == 'quadrangle':
                if fe_order == 'linear':
                    local_nodes[0,0,:] = [0,1,1,0]
                    local_nodes[0,1,:] = [0,0,1,1]
                elif fe_order == 'quadratic':
                    local_nodes[0,0,:] = [0,2,2,0,1,2,1,0,1]
                    local_nodes[0,1,:] = [0,0,2,2,0,1,2,1,1]
                elif fe_order == 'serendipity':
                    local_nodes[0,0,:] = [0,2,2,0,1,2,1,0]
                    local_nodes[0,1,:] = [0,0,2,2,0,1,2,1]
                else:
                    raise NotImplementedError(
                        f"Unsupported fe_order: {fe_order} for quadrangle in set_local_node_indexing!")

            elif fe_type == 'triangle':
                if fe_order == 'linear':
                    local_nodes[0,0,:] = [0,1,0]
                    local_nodes[0,1,:] = [0,0,1]
                    local_nodes[1,0,:] = [1,0,1]
                    local_nodes[1,1,:] = [1,1,0]
                elif fe_order == 'quadratic':
                    local_nodes[0,0,:] = [0,2,0,1,1,0]
                    local_nodes[0,1,:] = [0,0,2,0,1,1]
                    local_nodes[1,0,:] = [2,0,2,1,1,2]
                    local_nodes[1,1,:] = [2,2,0,2,1,1]
                else:
                    raise NotImplementedError(
                        f"Unsupported fe_order: {fe_order} for triangle in set_local_node_indexing!")

        elif ndim == 3:

            if fe_type == 'hexahedron':
                if fe_order == 'linear':
                    local_nodes[0,0,:] = [0,1,1,0,0,1,1,0]
                    local_nodes[0,1,:] = [0,0,1,1,0,0,1,1]
                    local_nodes[0,2,:] = [0,0,0,0,1,1,1,1]
                elif fe_order == 'quadratic':
                    local_nodes[0,0,:] = [0,2,2,0,0,2,2,0, 1,2,1,0,1,2,1,0, 0,2,2,0, 1,2,1,0, 1,1,1]
                    local_nodes[0,1,:] = [0,0,2,2,0,0,2,2, 0,1,2,1,0,1,2,1, 0,0,2,2, 0,1,2,1, 1,1,1]
                    local_nodes[0,2,:] = [0,0,0,0,2,2,2,2, 0,0,0,0,2,2,2,2, 1,1,1,1, 1,1,1,1, 0,2,1]
                elif fe_order == 'serendipity':
                    local_nodes[0,0,:] = [0,2,2,0,0,2,2,0, 1,2,1,0,1,2,1,0, 0,2,2,0]
                    local_nodes[0,1,:] = [0,0,2,2,0,0,2,2, 0,1,2,1,0,1,2,1, 0,0,2,2]
                    local_nodes[0,2,:] = [0,0,0,0,2,2,2,2, 0,0,0,0,2,2,2,2, 1,1,1,1]
                else:
                    raise NotImplementedError(
                        f"Unsupported fe_order: {fe_order} for hexahedron in set_local_node_indexing!")

            elif fe_type == 'tetrahedron':
                if fe_order == 'linear':
                    # Element 0: linear hexa nodes [0, 1, 2, 6]
                    local_nodes[0,0,:] = [0,1,1,1]
                    local_nodes[0,1,:] = [0,0,1,1]
                    local_nodes[0,2,:] = [0,0,0,1]
                    # Element 1: linear hexa nodes [0, 2, 3, 6]
                    local_nodes[1,0,:] = [0,1,0,1]
                    local_nodes[1,1,:] = [0,1,1,1]
                    local_nodes[1,2,:] = [0,0,0,1]
                    # Element 2: linear hexa nodes [0, 1, 6, 5]
                    local_nodes[2,0,:] = [0,1,1,1]
                    local_nodes[2,1,:] = [0,0,1,0]
                    local_nodes[2,2,:] = [0,0,1,1]
                    # Element 3: linear hexa nodes [0, 5, 6, 4]
                    local_nodes[3,0,:] = [0,1,1,0]
                    local_nodes[3,1,:] = [0,0,1,0]
                    local_nodes[3,2,:] = [0,1,1,1]
                    # Element 4: linear hexa nodes [0, 6, 3, 7]
                    local_nodes[4,0,:] = [0,1,0,0]
                    local_nodes[4,1,:] = [0,1,1,1]
                    local_nodes[4,2,:] = [0,1,0,1]
                    # Element 5: linear hexa nodes [0, 6, 7, 4]
                    local_nodes[5,0,:] = [0,1,0,0]
                    local_nodes[5,1,:] = [0,1,1,0]
                    local_nodes[5,2,:] = [0,1,1,1]

                elif fe_order == 'quadratic':
                    # 4 Corners (0,1,2,3), then 6 Mid-edges (0-1, 1-2, 0-2, 0-3, 1-3, 2-3)
                    # Element 0: linear hexa nodes [0, 1, 2, 6]
                    local_nodes[0,0,:] = [0,2,2,2, 1,2,1,1,2,2]
                    local_nodes[0,1,:] = [0,0,2,2, 0,1,1,1,1,2]
                    local_nodes[0,2,:] = [0,0,0,2, 0,0,0,1,1,1]
                    # Element 1: linear hexa nodes [0, 2, 3, 6]
                    local_nodes[1,0,:] = [0,2,0,2, 1,1,0,1,2,1]
                    local_nodes[1,1,:] = [0,2,2,2, 1,2,1,1,2,2]
                    local_nodes[1,2,:] = [0,0,0,2, 0,0,0,1,1,1]
                    # Element 2: linear hexa nodes [0, 1, 6, 5]
                    local_nodes[2,0,:] = [0,2,2,2, 1,2,1,1,2,2]
                    local_nodes[2,1,:] = [0,0,2,0, 0,1,1,0,0,1]
                    local_nodes[2,2,:] = [0,0,2,2, 0,1,1,1,1,2]
                    # Element 3: linear hexa nodes [0, 5, 6, 4]
                    local_nodes[3,0,:] = [0,2,2,0, 1,2,1,0,1,1]
                    local_nodes[3,1,:] = [0,0,2,0, 0,1,1,0,0,1]
                    local_nodes[3,2,:] = [0,2,2,2, 1,2,1,1,2,2]
                    # Element 4: linear hexa nodes [0, 6, 3, 7]
                    local_nodes[4,0,:] = [0,2,0,0, 1,1,0,0,1,0]
                    local_nodes[4,1,:] = [0,2,2,2, 1,2,1,1,2,2]
                    local_nodes[4,2,:] = [0,2,0,2, 1,1,0,1,2,1]
                    # Element 5: linear hexa nodes [0, 6, 7, 4]
                    local_nodes[5,0,:] = [0,2,0,0, 1,1,0,0,1,0]
                    local_nodes[5,1,:] = [0,2,2,0, 1,2,1,0,1,1]
                    local_nodes[5,2,:] = [0,2,2,2, 1,2,1,1,2,2]
                else:
                    raise NotImplementedError(
                        f"Unsupported fe_order: {fe_order} for tetrahedron in set_local_node_indexing!")

        else:
            raise NotImplementedError(f"Unsupported ndim ={ndim} in set_local_node_indexing!")

        return local_nodes

    def set_rectangular_mesh_serendipity(self, nel_dim: NDArray[np.int32],
        nel_per_cell: int, major_order: str = 'row'):
        """Sets up a rectangular structured FE mesh for serendipity elements.

        Args:
            self: the MeshDiscretization object that holds (nbf, fe_type)
            nel_dim: the array of the number of elements in each dimension
            nel_per_cell: the number of elements in which a cell is divided
            major_order (optional): the order of enumeration
        Returns:
            The structure FE mesh data in self (nel, x_mesh, connectivity,
            boundary_nodes)."""

        fe_type = self.fe_type
        nbf0 = self.nbf
        nel = self.nel

        # 'quadratic' to generate regular connectivity
        if fe_type == 'quadrangle':
            nbf = 9
            center_idx = [8]
        elif fe_type == 'hexahedron':
            nbf = 27
            center_idx = np.array([20, 21, 22, 23, 24, 25, 26], dtype=int)
        else:
            raise NotImplementedError('Unsupported fe_type in set_rectangular_mesh_serendipity!')

        self.fe_order = 'quadratic'
        self.nbf = nbf
        self.set_rectangular_mesh_connectivity_fe(nel_dim, nel_per_cell, major_order)

        # center nodes global IDs
        center_nodes = set()
        for iel in range(nel):
            for i in range(len(center_idx)):
                idx = center_idx[i]
                center_nodes.add(self.connectivity[iel, idx])

        # Filter the mesh coordinate array and create an index map
        x_new = []
        old_to_new_mapping = {}
        new_idx = 0
        for old_idx in range(self.nodes_total):
            if old_idx in center_nodes:
                continue
            else:
                x_new.append(self.x_mesh[old_idx])
                old_to_new_mapping[old_idx] = new_idx
                new_idx += 1

        self.x_mesh = np.array(x_new)
        self.nodes_total = self.x_mesh.shape[0]

        # Create the serendipity connectivity
        self.fe_order = 'serendipity'
        self.nbf = nbf0
        self.set_rectangular_mesh_connectivity_fe(nel_dim, nel_per_cell, major_order)

        for iel in range(nel):
            for inod in range(self.nbf):
                old_gnod = int(self.connectivity[iel, inod])
                self.connectivity[iel, inod] = old_to_new_mapping[old_gnod]

    def set_rectangular_mesh_node_types_fe(self):
        """Classifies the nodes for a 2D FE rectangular mesh (boundary, multi-boundary).

        Args:
            self: the MeshDiscretization object that holds (nodes_total)
        Returns:
            The lists of node types in self for each type of node
            (boundary_nodes)."""

        nodes_total = self.nodes_total
        ndim = self.ndim
        x0 = self.x0
        xf = self.xf
        x_mesh = self.x_mesh

        self.bulk_nodes = []
        self.boundary_nodes = []

        for i in range(nodes_total):

            xi = x_mesh[i,:]

            boundaries = []
            for d in range(ndim):
                if np.abs(xi[d] - x0[d]) < 1e-8:
                    bnd = 2*d
                    boundaries.append(bnd)
                if np.abs(xi[d] - xf[d]) < 1e-8:
                    bnd = 2*d+1
                    boundaries.append(bnd)

            bnd_count = len(boundaries)
            if bnd_count == 0:
                self.bulk_nodes.append(i)
            else:
                self.boundary_nodes.append((i, boundaries))
