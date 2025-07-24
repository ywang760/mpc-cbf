import numpy as np


def generate_points_on_circle(num_points, radius, angle_bias):
    """
    Generates (num_points) points uniformly at random on the circumference
    of a circle with a given (radius).

    Returns:
        x_coords (ndarray): Array of x-coordinates.
        y_coords (ndarray): Array of y-coordinates.
    """
    # Incremental angles from [angle_bias, angle_bias+2pi]
    thetas = np.linspace(angle_bias, angle_bias+2*np.pi, num=num_points, endpoint=False)

    # Convert polar to Cartesian
    x_coords = radius * np.cos(thetas)
    y_coords = radius * np.sin(thetas)

    # make sure the points are rounded to 3 decimal places
    x_coords = np.round(x_coords, 3)
    y_coords = np.round(y_coords, 3)

    return x_coords, y_coords

def compute_yaw(x, y, circle_center):
    yaws = np.arctan2(-y - circle_center[1], -x - circle_center[0])
    return yaws

def generate_points_on_formation(num_points, n_row, distance_x, distance_y, start_bias):
    n_col = np.ceil(num_points/n_row)

    x_min = -(n_col/2)*distance_x
    x_max = (n_col/2)*distance_x
    y_min = -(n_row/2)*distance_y
    y_max = (n_row/2)*distance_y

    x = np.linspace(x_max, x_min, int(n_col))
    y = np.linspace(y_max, y_min, int(n_row))

    xx, yy = np.meshgrid(x, y)

    # all_points_in_formation = np.hstack([xx.flatten().reshape(-1,1), yy.flatten().reshape(-1,1)])
    x_coords = []
    y_coords = []
    for i in range(num_points):
        x_coords.append(xx.flatten()[i] + start_bias[0])
        y_coords.append(yy.flatten()[i] + start_bias[1])
    x_coords = np.array(x_coords)
    y_coords = np.array(y_coords)
    return x_coords, y_coords


def plot_position(x, y, idx, c="r", ax=None):
    ax.scatter(x, y, c=c)
    ax.text(x, y, str(idx), fontsize=10)
