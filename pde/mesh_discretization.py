
import numpy as np

class MeshDiscretization():

    def __init__(self, x0: np.ndarray[tuple[int]], xf: np.ndarray[tuple[int]],
                nodes_dim: np.ndarray[tuple[int]] = None,
                x_mesh: np.ndarray[tuple[int, int]] = None,
                connectivity: np.ndarray[tuple[int, int]] = None,
                boundary_nodes: np.ndarray[tuple[int]] = None):
        
        self.x0 = np.array(x0, dtype=float)
        self.xf = np.array(xf, dtype=float)

        ndim = x0.shape[0]
        self.ndim = ndim

        if nodes_dim is None:
            self.nodes_dim = np.full(ndim, 5, dtype=int)
        else:
            self.nodes_dim = np.array(nodes_dim, dtype=int)

        self.nodes_total = 1
        for d in range(ndim):
            self.nodes_total *= self.nodes_dim[d]

        if x_mesh is None:
            self.x_mesh = np.zeros((self.nodes_total,ndim))
        else: 
            self.x_mesh = x_mesh
        
        self.connectivity = connectivity

        self.boundary_nodes = boundary_nodes

    def read_mesh(self):
        pass

    def read_connectivity(self):
        pass
    
    def set_rectangular_mesh(self, p: np.ndarray[tuple[int]] = None,
        major_order: str = 'row'):

        self.x_mesh = self.create_rectangular_mesh(p, major_order)
        self.connectivity = self.structured_connectivity(major_order)
        self.bulk_nodes, self.boundary_nodes, self.multi_boundary_nodes = self.classify_nodes()

    def create_rectangular_mesh(self, p: np.ndarray[tuple[int]] = None,
        major_order: str = 'row'):
        """Sets the mesh coordinates for a rectangular mesh
        given the stretching factor for refinement in each direction
        and the order of node enumeration.

        Args:
            self: the MeshDiscretization object that holds (x0, xf, ndim)
            p (optional): the power law refinement factor.
                For p > 1: refinement towards x0.
                For 0 < p < 1: refinement towards xf
            major_order (optional): the order of enumeration
        Returns:
            The arrays of mesh coordinates."""

        ndim = self.ndim

        if p is None:
            p = np.ones(ndim, dtype=float)

        x_1d = []
        for d in range(ndim):
            xd = np.array([self.x0[d], self.xf[d]])
            nd = self.nodes_dim[d]
            x_1d.append(self.create_mesh_1d(nd, xd, p[d]))

        if major_order == 'row':
            # row-major order: the last dimension (z) cycles first
            dim_start, dim_end, dim_step = ndim-1, -1, -1
        else:
            # column-major order: the first dimension (x) cycles first
            dim_start, dim_end, dim_step = 0, ndim, +1
        
        x_mesh = np.zeros((self.nodes_total,ndim))
        for inod in range(self.nodes_total):
            idx = inod
            for d in range(dim_start, dim_end, dim_step):
                nd = self.nodes_dim[d]
                idx_dim = idx%nd
                x_mesh[inod,d] = x_1d[d][idx_dim]
                idx = idx//nd

        return x_mesh

    def create_mesh_1d(self, nodes_dim: int, xd: np.ndarray[tuple[int]],
        p: int = 1) -> np.ndarray[tuple[int]]:
        """Returns the mesh coordinates given the number of nodes,
        the domain bounds and the stretching factor for refinement.

        Args:
            nodes_dim: the number of nodes
            xd: the domain bounds
            p (optional): the power law refinement factor.
                For p > 1: refinement towards x0.
                For 0 < p < 1: refinement towards xf
        Returns:
            The array of 1D mesh coordinates."""

        x0, xf = xd[0], xd[1]
        L = xf - x0

        x = np.zeros(nodes_dim)
        for i in range(nodes_dim):
            xi = i/(nodes_dim-1)
            gi = xi**p
            # gi = 0.5*(1.0-np.cos(p*np.pi*xi))
            # gi = np.tanh(p*(2*xi-1))/(2*np.tanh(p)) + 0.5
            x[i] = x0 + L*gi

        return x
    
    def structured_connectivity(self, major_order: str):
        """Sets the connectivity for a rectangular mesh
        given the order of node enumeration.

        Args:
            self: the MeshDiscretization object that holds (ndim, nodes_dim, nodes_total)
            major_order (optional): the order of enumeration
        Returns:
            The connectivity array."""

        ndim = self.ndim
        nodes_dim = self.nodes_dim
        nodes_total = self.nodes_total

        node_step = [1]*ndim
        if major_order == 'row':
            # row-major order: the last dimension (z) cycles the fastest
            for d in range(ndim - 2, -1, -1):
                node_step[d] = node_step[d+1]*nodes_dim[d+1]
            dim_start, dim_end, dim_step = ndim-1, -1, -1
        else:
            # column-major order: the first dimension (x) cycles the fastest
            for d in range(1, ndim):
                node_step[d] = node_step[d-1]*nodes_dim[d-1]
            dim_start, dim_end, dim_step = 0, ndim, +1
        
        connectivity = np.full((nodes_total, ndim, 2), -1, dtype=int)

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
                    connectivity[inod, d, 0] = inod - node_step[d]
                    
                # Forward Neighbor check
                if idx_dim[d] < nodes_dim[d] - 1:
                    connectivity[inod, d, 1] = inod + node_step[d]

        return connectivity
    
    def classify_nodes(self):
        """Classifies the nodes for a rectangular mesh (bulk, boundary, multi-boundary).

        Args:
            self: the MeshDiscretization object that holds (nodes_total, ndim, connectivity)
        Returns:
            A tuple of arrays for each type of node."""

        nodes_total = self.nodes_total
        ndim = self.ndim
        connectivity = self.connectivity

        bulk_nodes = []
        boundary_nodes = []
        # boundary_nodes = np.full((nodes_total,ndim), fill_value=-1, dtype=int)
        multi_boundary_nodes = []
        
        for i in range(nodes_total):

            bnd_count = 0
            boundaries = []
            # the highest dimension index dominates corners
            for d in range(ndim):
                for j in range(connectivity.shape[2]):
                    if connectivity[i, d, j] == -1:
                        bnd_count += 1
                        bnd = 2*d+j
                        boundaries.append(bnd)
                        break

            if bnd_count == 0:
                bulk_nodes.append(i)
            elif bnd_count == 1:
                boundary_nodes.append((i, boundaries))
            else:
                multi_boundary_nodes.append((i, boundaries))

        return bulk_nodes, boundary_nodes, multi_boundary_nodes
    
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

        # dx = np.diff(self.x)
        dx = self.x[1:] - self.x[0:-1]
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

        # dx = np.diff(self.x)
        dx = self.x[1:] - self.x[0:-1]
        dx_max = np.max(dx)

        pe = (v*dx_max)/D

        dx_stable = D*2.0/v
        if pe > 2.0:
            print(f"dx = {dx_max:.4f} > dx_max = {dx_stable:.4f} (Pe = {pe:.1f} > 2.0)")
            raise ValueError('Too large dx for this system')
