# -*- coding: utf-8 -*-
import numpy as np

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


# def octahedral_coordinate(nt: int):

#     # Do the (x + y + z = nt) face of the octahedron
#     # z -> 0 to nt-1
#     # y -> 0 to nt-z
#     # x -> nt - y - z

#     n = int((nt + 1) * (nt + 2) / 2)
#     x = np.empty(n, dtype=np.float64)
#     y = np.empty(n, dtype=np.float64)
#     z = np.empty(n, dtype=np.float64)

#     k = 0
#     for j in range(nt + 1):
#         for i in range(nt - j + 1):
#             # x = nt-i-j;
#             # y = i;
#             # z = j;
#             x[k] = nt - i - j
#             y[k] = i
#             z[k] = j
#             k += 1

#     return x, y, z


def octahedral_direction_cosine_squares_and_amplitudes(nt: int):

    # Do the (x + y + z = nt) face of the octahedron
    # z -> 0 to nt-1
    # y -> 0 to nt-z
    # x -> nt - y - z

    n = int((nt + 1) * (nt + 2) / 2)
    x = np.empty(n, dtype=np.float64)
    y = np.empty(n, dtype=np.float64)
    z = np.empty(n, dtype=np.float64)

    k = 0
    for j in range(nt + 1):
        for i in range(nt - j + 1):
            # x = nt-i-j;
            # y = i;
            # z = j;
            x[k] = nt - i - j
            y[k] = i
            z[k] = j
            k += 1

    x *= x
    y *= y
    z *= z

    r2 = x + y + z

    x /= r2
    y /= r2
    z /= r2
    amp = nt / (r2 * np.sqrt(r2))

    return x, y, z, amp


def cosine_of_polar_angles_and_amplitudes(integration_density: int = 72):
    r"""
    Calculates and return the direction cosines and the related amplitudes for
    the positive quadrant of the sphere. The direction cosines corresponds to
    angle $\alpha$ and $\beta$, where $\alpha$ is the azimuthal angle and
    $\beta$ is the polar angle. The amplitudes are evaluated as $\frac{1}{r^3}$
    where $r$ is the distance from the origin to the face of the unit
    octahedron in the positive quadrant along the line given by the values of
    $\alpha$ and $\beta$.

    :ivar integration_density:
        The value is an integer which represents the frequency of class I
        geodesic polyhedra. These polyhedra are used in calculating the
        spherical average. Presently we only use octahedral as the frequency1
        polyhedra. As the frequency of the geodesic polyhedron increases, the
        polyhedra approach a sphere geometry. A higher frequency will result in
        a better powder averaging. The default value is 72.
        Read more on the `Geodesic polyhedron
        <https://en.wikipedia.org/wiki/Geodesic_polyhedron>`_.

    :return cos_alpha: The cosine of the azimuthal angle.
    :return cos_beta: The cosine of the polar angle.
    :return amp: The amplitude at the given $\alpha$ and $\beta$.
    """
    nt = integration_density
    xr, yr, zr, amp = octahedral_direction_cosine_squares_and_amplitudes(nt)

    cos_beta = np.sqrt(zr)
    cos_alpha = np.zeros(xr.size, dtype=np.float64)
    cos_alpha[:-1] = np.sqrt(xr[:-1] / (xr[:-1] + yr[:-1]))
    cos_alpha[-1] = 0.0
    return cos_alpha, cos_beta, amp


def right_triangle_process(spec, l_clip_rt, r_clip_rt, r_clip_lt, top, f21, pmax, p, f):
    if p != pmax:
        p = bin_right_half(spec, l_clip_rt, r_clip_rt, top, f21, pmax, p, f)
    else:
        if not r_clip_rt:
            spec[p] += f21 * top * 0.5
        if r_clip_rt and not r_clip_lt:
            df2 = top / f21
            diff = f[2] - p - 1.0
            spec[p] += (f21 - diff) * (diff + f21) * 0.5 * df2
    return p


def left_triangle_process(spec, l_clip_lt, r_clip_lt, l_clip_rt, top, f10, pmid, p, f):
    if p != pmid:
        p = bin_left_half(spec, l_clip_lt, r_clip_lt, top, f10, pmid, p, f)
    else:
        if not r_clip_lt and not l_clip_lt:
            spec[p] += f10 * top * 0.5
        # if l_clip_lt: # and not l_clip_rt:
        #     spec[p] += f[1] * (f10 - f[0]) * 0.5 * top / f10
    return p


def triangle_interpolation1D(f, spec, amp=1.0):
    """Interpolate between three points to form a triangle."""

    points = spec.size
    f = np.asarray(f, dtype=np.float64)

    p = int(f[0])
    if int(f[0]) == int(f[1]) == int(f[2]):
        if p >= points or p < 0:
            return
        spec[p] += amp
        return

    f = np.sort(f)

    top = amp * 2.0 / (f[2] - f[0])

    p = int(f[0])
    pmid = int(f[1])
    pmax = int(f[2])

    f10 = f[1] - f[0]
    f21 = f[2] - f[1]

    if pmax < 0 or p >= points:
        return

    clips, p, pmid, pmax = get_clip_conditions(f, p, pmid, pmax, points)
    l_clip_lt, l_clip_rt, r_clip_lt, r_clip_rt = clips
    print(clips)

    if f[1] >= 0:
        p = left_triangle_process(
            spec, l_clip_lt, r_clip_lt, l_clip_rt, top, f10, pmid, p, f
        )
    if f[2] >= 0:
        p = right_triangle_process(
            spec, l_clip_rt, r_clip_rt, r_clip_lt, top, f21, pmax, p, f
        )


def get_clip_conditions(f, p, pmid, pmax, points):
    """Return True for every clip boundary"""
    r_clip_lt = r_clip_rt = False
    l_clip_lt = l_clip_rt = False

    if p <= 0:
        p = 0
        l_clip_lt = True

    if pmid > points - 1:
        pmid = points - 1
        r_clip_lt = True

    if pmid <= 0:
        pmid = 0
        l_clip_rt = True

    if pmax > points - 1:
        pmax = points - 1
        r_clip_rt = True

    return [l_clip_lt, l_clip_rt, r_clip_lt, r_clip_rt], p, pmid, pmax


def bin_left_half(spec, l_clip_lt, r_clip_lt, top, f10, pmid, p, f):
    """Interpolate the left half of the triangle"""
    df1 = top / f10
    diff = p + 1.0 - f[0]

    spec[p] += (
        (diff - 0.5) * df1 if l_clip_lt and f[0] <= 0.0 else 0.5 * diff * diff * df1
    )

    p += 1

    diff -= 0.5
    diff *= df1
    while p != pmid:
        diff += df1
        spec[p] += diff
        p += 1

    p = pmid
    if not r_clip_lt:
        spec[p] += (f[1] - p) * (f10 + p - f[0]) * 0.5 * df1
    else:
        spec[p] += diff + df1
    return p


def bin_right_half(spec, l_clip_rt, r_clip_rt, top, f21, pmax, p, f):
    """Interpolate the right half of the triangle"""
    df1 = top / f21
    diff = f[2] - p - 1.0

    spec[p] += (
        (f21 - diff) * (diff + f21) * 0.5 * df1 if not l_clip_rt else (diff + 0.5) * df1
    )

    p += 1

    diff += 0.5
    diff *= df1
    while p != pmax:
        diff -= df1
        spec[p] += diff
        p += 1

    p = pmax
    if not r_clip_rt:
        spec[p] += (f[2] - p) ** 2 * 0.5 * df1
    else:
        spec[p] += diff - df1
    return p


# def average_over_octant(spec, freq, nt, amp):

#     n_pts = (nt + 1) * (nt + 2) / 2

#     #   Interpolate between frequencies by setting up tents

#     local_index = nt - 1
#     i = 0
#     j = 0
#     while i < n_pts - 1:
#         temp = amp[i + 1] + amp[nt + 1 + j]
#         amp1 = temp
#         amp1 += amp[i]

#         freq_3 = [freq[i], freq[i + 1], freq[nt + 1 + j]]
#         triangle_interpolation1D(freq_3, spec, amp1)

#         if i < local_index:
#             amp1 = temp
#             amp1 += amp[nt + 1 + j + 1]
#             freq_3 = [freq[i + 1], freq[nt + 1 + j], freq[nt + 1 + j + 1]]
#             triangle_interpolation1D(freq_3, spec, amp1)
#         else:
#             local_index = j + nt
#             i += 1
#         i += 1
#         j += 1
