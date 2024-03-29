import numpy as np
from mrsimulator.simulator.config import CustomSampling
from scipy.spatial import ConvexHull


def generate_custom_sampling(alpha, beta, weight, triangle_mesh=False):
    if not triangle_mesh:
        return CustomSampling(
            alpha=alpha, beta=beta, weight=weight, vertex_indexes=None
        )

    z = np.cos(beta)
    x = np.sin(beta) * np.cos(alpha)
    y = np.sin(beta) * np.sin(alpha)

    points = np.array([x, y, z])
    hull = ConvexHull(points.T)
    vertex_indexes = hull.simplices

    sampling = CustomSampling(
        alpha=alpha, beta=beta, weight=weight, vertex_indexes=vertex_indexes
    )
    return sampling


def step_averaging(N_alpha: int, N_beta: int, triangle_mesh=True):
    """Generate STEP averaging samples.

    Args:
        N_alpha: number of points along the alpha dimension.
        beta: number of points along the beta dimension.
        triangle_mesh: generate ConvexHull triangulation of points.
    """
    sphereType = "sphere"
    norm_step = 0.0
    if sphereType == "sphere":
        a = [1.0, 0.0, 1.0]
        b = 1.0
    if sphereType == "hemisphere":
        a = [1.0, 0.0, 1.0]
        b = 2.0
    if sphereType == "octant":
        a = [2.0, 1.0, 8.0]
        b = 2.0

    inc_alpha = 2 * np.pi / ((N_alpha - 1) * a[2])
    inc_beta = np.pi / (2.0 * (N_beta) * b)

    norm_step = np.sin(inc_beta * (2.0 * np.arange(N_beta) + 1.0)).sum()
    norm_step = 1.0 / (N_alpha * norm_step)

    alx, btx = np.meshgrid(np.arange(N_alpha), np.arange(N_beta), indexing="xy")
    alx = alx.ravel()
    btx = btx.ravel()

    alpha = inc_alpha * (a[1] + alx * a[0])
    beta = inc_beta * (2.0 * btx + 1.0)
    weight = norm_step * np.sin(inc_beta * (2.0 * btx + 1.0))
    # shift alpha by pi/2 for beta > pi/2
    alpha[beta > np.pi / 2] += (np.pi / 2) / (2 * N_alpha)
    return generate_custom_sampling(alpha, beta, weight, triangle_mesh)


def get_zcw_number(M):
    """ZCW number"""
    # returns the number of ZCW angles for the given integer M=2,3,4,...
    gM = 5
    gMminus1 = 3
    sum = 5
    for _ in range(M + 1):
        sum = gM + gMminus1
        gMminus1 = gM
        gM = sum
    return sum


def zcw_averaging(M: int, triangle_mesh=True):
    """Generate ZCW averaging samples.

    Args:
        M: ZCW point generation factor.
        triangle_mesh: generate ConvexHull triangulation of points.
    """
    sphereType = "sphere"
    if sphereType == "sphere":
        c = [1.0, 2.0, 1.0]
    if sphereType == "hemisphere":
        c = [-1.0, 1.0, 1.0]
    if sphereType == "octant":
        c = [-1.0, 1.0, 4.0]

    N = get_zcw_number(M)
    g2 = get_zcw_number(M - 2)

    m = np.arange(N)
    beta = np.arccos(c[0] * (c[1] * np.fmod(m / float(N), 1.0) - 1.0))
    alpha = 2.0 * np.pi * (np.fmod(m * float(g2) / float(N), 1.0)) / c[2]
    weight = np.ones(N, dtype=float) / N
    return generate_custom_sampling(alpha, beta, weight, triangle_mesh)


if __name__ == "__main__":
    sampling = zcw_averaging(M=21)
    rad2deg = 180.0 / np.pi
    nd_array = np.array(
        [sampling.alpha * rad2deg, sampling.beta * rad2deg, sampling.weight]
    ).T
    size = sampling.alpha.size
    np.savetxt(f"zcw{size}.cry", nd_array, header=str(size), fmt="%.6e")
