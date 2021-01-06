# -*- coding: utf-8 -*-
"""
Example Script:  Nonconvexity of Gaussians Demo
"""

import matplotlib.pyplot as plt
import matplotlib.cm
import scipy
import scipy.stats
import numpy as np


def get_grid_points(lb, ub, num_pts):
    return lb + (ub - lb) * (np.arange(num_pts, dtype=float)) / (num_pts - 1)


def get_gaussians(arr, mean, stdev):
    """
    returns an array of gaussian pdf values, which have been normalized so that
    the maximum of the pdf is one.
    """

    z = (arr - mean) / stdev
    #    print(arr,z)
    arr = scipy.stats.norm.pdf(z) / stdev
    return arr


def get_area_included(lb, ub, mean, stdev):
    zl = (lb - mean) / stdev
    zu = (ub - mean) / stdev
    return scipy.stats.norm.cdf(zu) - scipy.stats.norm.cdf(zl)


def feasibility_check(x, bounds):
    return (x > bounds).any()


def plot_results_vector(x, vectors, names, fname):
    plt.cla()
    plt.clf()
    f, ax = plt.subplots(nrows=1, ncols=1)
    for v in range(len(vectors)):
        # Question from Alex M: Did you mean to do 'ax.plt.plot(...) here?
        plt.plot(list(x), list(vectors[v]), label=names[v])
    plt.legend()
    plt.savefig(fname)


#    plt.show()

def search_best_pair(bounds, s1, s2, pts_per_dim, est_points, lb=0., ub=1.):
    best_obj = 0.
    best_sol = (0., 0.)
    x_search = np.arange(pts_per_dim, dtype=float) / (pts_per_dim - 1)
    #    print(x1_search,x2_search)
    objs = np.zeros([pts_per_dim, pts_per_dim], dtype=float)
    est_objs = np.zeros([pts_per_dim, pts_per_dim], dtype=float)
    grid = get_grid_points(lb, ub, len(bounds))
    for i in range(len(x_search)):
        for j in range(len(x_search)):
            x = get_gaussians(grid, x_search[i], s1) + get_gaussians(grid, x_search[j], s2)
            #            print(x1_search[i],x2_search[j],x,grid,lb,ub)
            if not (x > bounds).any():  # feasibility check
                objs[i, j] = get_area_included(lb, ub, x_search[i], s1) + get_area_included(lb, ub, x_search[j], s2)
                est_objs[i, j] += (
                                          get_gaussians(est_points, x_search[i], s1) +
                                          get_gaussians(est_points, x_search[j], s2)
                                  ).sum() / len(est_points)
                if objs[i, j] > best_obj:
                    best_obj = objs[i, j]
                    best_sol = (x_search[i], x_search[j])
                #                assert(False)
    return best_sol, best_obj, objs, est_objs


def search_best_triple(bounds, s1, s2, s3, pts_per_dim, est_points, lb=0., ub=1.):
    best_obj = 0.
    best_sol = (0., 0., 0.)
    x_search = np.arange(pts_per_dim, dtype=float) / (pts_per_dim - 1)
    #    print(x1_search,x2_search)
    objs = np.zeros([pts_per_dim, pts_per_dim, pts_per_dim], dtype=float)
    est_objs = np.zeros([pts_per_dim, pts_per_dim, pts_per_dim], dtype=float)
    grid = get_grid_points(lb, ub, len(bounds))
    for i in range(len(x_search)):
        for j in range(len(x_search)):
            for k in range(len(x_search)):
                x = get_gaussians(grid, x_search[i], s1) + get_gaussians(grid, x_search[j], s2) + get_gaussians(grid,
                                                                                                                x_search[
                                                                                                                    k],
                                                                                                                s3)
                #            print(x1_search[i],x2_search[j],x,grid,lb,ub)
                if not (x > bounds).any():  # feasibility check
                    objs[i, j, k] = get_area_included(lb, ub, x_search[i], s1) + get_area_included(lb, ub, x_search[j],
                                                                                                   s2) + get_area_included(
                        lb, ub, x_search[k], s3)
                    est_objs[i, j, k] += (
                                                 get_gaussians(est_points, x_search[i], s1) +
                                                 get_gaussians(est_points, x_search[j], s2) +
                                                 get_gaussians(est_points, x_search[k], s3)
                                         ).sum() / len(est_points)
                    if objs[i, j, k] > best_obj:
                        best_obj = objs[i, j, k]
                        best_sol = (x_search[i], x_search[j], x_search[k])
                        #                assert(False)
    return best_sol, best_obj, objs, est_objs


def plot_obj_heatmap(objs, fname):
    plt.cla()
    plt.clf()
    plt.imshow(objs, cmap='hot')
    plt.colorbar()
    plt.savefig(fname)


def plot_2d_results(best_sol, best_obj, objs, est_objs, bounds, lb, ub, num_pts, fname, s1, s2):
    """
    Function: create a plot of each of the results, along with a legend
    :param best_sol:
    :param best_obj:
    :param objs:
    :param est_objs:
    :param bounds:
    :param lb: Lower bound
    :param ub: Upper bound
    :param num_pts:
    :param fname: output file name (as a pdf)
    :param s1:
    :param s2:
    :return:
    """
    x = get_grid_points(lb, ub, num_pts)
    v1 = get_gaussians(x, best_sol[0], s1)
    v2 = get_gaussians(x, best_sol[1], s2)
    v3 = v1 + v2
    plot_results_vector(x, [v1, v2, v3, bounds], ["x1", "x2", "x1+x2", "bounds"], fname)


def plot_3d_results(best_sol, best_obj, objs, est_objs, bounds, lb, ub, num_pts, fname, s1, s2, s3):
    """
    Function: Plots multiple gaussians, along with the sum of them
    :param best_sol:
    :param best_obj:
    :param objs:
    :param est_objs:
    :param bounds:
    :param lb: Lower bound
    :param ub:
    :param num_pts:
    :param fname: output file name (pdf)
    :param s1:
    :param s2:
    :param s3:
    :return:
    """
    x = get_grid_points(lb, ub, num_pts)
    v1 = get_gaussians(x, best_sol[0], s1)
    v2 = get_gaussians(x, best_sol[1], s2)
    v3 = get_gaussians(x, best_sol[2], s3)
    v4 = v1 + v2 + v3
    plot_results_vector(x, [v1, v2, v3, v4, bounds], ["x1", "x2", "x3", "x1+x2+x3", "bounds"], fname)


def plot_optimal_flux_heatmap(outputs,fname):
    flux = outputs.flux_map.reshape([outputs.flux_model.receiver.params["pts_per_dim"],outputs.flux_model.receiver.params["pts_per_dim"]])
    plt.imshow(flux, cmap='hot')
    plt.colorbar()
    plt.savefig(fname)
    plt.cla()
    plt.clf()
    
def plot_flux_violation(outputs,fname):
    """
    Plots heatmap of flux_violation at each measurement point. A value of zero is assigned to those
    points at which there is no flux-violation 

    Parameters
    ----------
    outputs : output - results - from the optimization model - Contains results along with objects of flux,field, receiver etc
    fname : TYPE
        DESCRIPTION.

    Returns
    -------
    Plots heatmap

    """
    pts = outputs.flux_model.receiver.params["pts_per_dim"]
    flux_ub = outputs.flux_model.receiver.flux_upper_limits
    flux = outputs.flux_map
    flux_violation = np.zeros_like(flux)
    for m in range(len(flux)):
        flux_violation[m] = max(0.0, flux[m]-flux_ub[m])
    flux_violation = np.array(flux_violation).reshape(pts,pts)
    plt.imshow(flux_violation, cmap='hot')
    plt.colorbar()
    plt.savefig(fname)
    plt.cla()
    plt.clf()
    ## The following code plots column wise max-flux-violation. 
    # max_col = np.max(flux_violation, axis = 0)
    # max_col = max_col.reshape(1,pts)
    # plt.imshow(max_col, cmap='hot')
    # plt.colorbar()
    # plt.savefig(fname+"_col_max")
    # plt.cla()
    # plt.clf()
    
def plot_field(outputs,fname):
    x = outputs.flux_model.field.x.flatten()
    y = outputs.flux_model.field.y.flatten()
    plt.scatter(x, y, s=6)
    plt.tight_layout()
    plt.savefig(fname, dpi = 2000)
    
def plot_optimal_aimpoint_allocation(outputs,fname):
    x = outputs.flux_model.field.x.flatten()
    y = outputs.flux_model.field.y.flatten()
    colors = plt.scatter(x,y,s=6,c=outputs.aimpoint_select_map,cmap='hsv')
    plt.savefig(fname, dpi= 2000)
    plt.cla()
    plt.clf()
    
def plot_optimal_aimpoint_guide(outputs,fname):
    x = outputs.flux_model.receiver.aim_x.flatten()
    y = outputs.flux_model.receiver.aim_z.flatten()
    colors = plt.scatter(x,y,s=6,c=range(1,outputs.flux_model.receiver.num_aimpoints+1),cmap='hsv')
    plt.savefig(fname, dpi= 2000)
    plt.cla()
    plt.clf()

def plot_defocused(outputs,fname):
    x = outputs.flux_model.field.x.flatten()
    y = outputs.flux_model.field.y.flatten()
    color = []
    for i in range(len(outputs.aimpoint_select_map)):
        color.append('black' if outputs.aimpoint_select_map[i] == 0 else 'r')
    colors = plt.scatter(x,y,s=5,c=color,cmap='hsv')
    plt.savefig(fname, dpi = 2000)
    plt.cla()
    plt.clf()
    
def plot_field_sections(outputs, fname):
    x = []
    y = []
    col = []    #Color
    for idx in range(len(outputs.flux_model.field.helios_by_section)):
            #r = np.random.random_sample()
        for h in outputs.flux_model.field.helios_by_section[idx]:
            x.append(outputs.flux_model.field.x.flatten()[h])
            y.append(outputs.flux_model.field.y.flatten()[h])
            col.append(float(idx))
    plt.scatter(x,y,s=5,c = col,cmap="jet")
    plt.savefig(fname, dpi= 2000)
    plt.cla()
    plt.clf()
if __name__ == "__main__":
    lb = 0.0
    ub = 1.0
    s1 = 0.1
    s2 = 0.2
    s3 = 0.4
    pts_per_dim = 51
    num_pts = 1001
    est_points = get_grid_points(lb, ub, 11)[:-1] + 0.05

    bounds = np.ones(num_pts, dtype=float) * 4. + 1. * np.arange(num_pts, dtype=float) / (num_pts - 1)
    best_sol, best_obj, objs, est_objs = search_best_pair(bounds, s1, s2, pts_per_dim, est_points, lb, ub)
    #    plot_2d_results(best_sol,best_obj,objs,est_objs,bounds,lb,ub,num_pts,"2d_results.pdf",s1,s2)
    print("2d results complete and plotted, sol:", best_sol)
    best_sol, best_obj, objs, est_objs = search_best_triple(bounds, s1, s2, s3, pts_per_dim, est_points, lb, ub)
    print("3d results complete, sol:", best_sol)
    plot_3d_results(best_sol, best_obj, objs, est_objs, bounds, lb, ub, num_pts, "3d_results_case2.pdf", s1, s2, s3)
    print("3d results plotted, sol:", best_sol)

    bounds = np.ones(num_pts, dtype=float) * 5. + 1.5 * np.arange(num_pts, dtype=float) / (num_pts - 1)
    best_sol, best_obj, objs, est_objs = search_best_triple(bounds, s1, s2, s3, pts_per_dim, est_points, lb, ub)
    print("3d results case 3 complete, sol:", best_sol)
    plot_3d_results(best_sol, best_obj, objs, est_objs, bounds, lb, ub, num_pts, "3d_results_case3.pdf", s1, s2, s3)
    bounds = np.ones(num_pts, dtype=float) * 6 + 2. * np.arange(num_pts, dtype=float) / (num_pts - 1)
    best_sol, best_obj, objs, est_objs = search_best_triple(bounds, s1, s2, s3, pts_per_dim, est_points, lb, ub)
    print("3d results case 4 complete, sol:", best_sol)
    plot_3d_results(best_sol, best_obj, objs, est_objs, bounds, lb, ub, num_pts, "3d_results_case4.pdf", s1, s2, s3)
    bounds = np.ones(num_pts, dtype=float) * 10 + 0 * np.arange(num_pts, dtype=float) / (num_pts - 1)
    best_sol, best_obj, objs, est_objs = search_best_triple(bounds, s1, s2, s3, pts_per_dim, est_points, lb, ub)
    print("3d results case 5 complete, sol:", best_sol)
    plot_3d_results(best_sol, best_obj, objs, est_objs, bounds, lb, ub, num_pts, "3d_results_case5.pdf", s1, s2, s3)

#    print(v1)
#    print(v2)
#    print(x)