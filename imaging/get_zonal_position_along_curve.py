import sys
import cv2
import numpy as np
import json
from functools import partial
import pandas as pd
import matplotlib.pyplot as plt

# Script used to annotate the position of GFP+ cells (teto-Olfr17 expressing cells) relative to a curve drawn on the most basal side of the MOE
# the curve passes through each of the 5 zones (starting in zone1 and ending in zone5) allowing to quantify expression in each zone
# used for Figure 7E

# input is immunofluorescnece images of GFP expression in the MOE
# script is run in stages:
# 1. in 'line' mode for annotation of the curve
# 2. in 'cell' model for annoation of the GFP+ cells
# 3. in 'analyze' mode, for combining the results of the earlier stages and producing a TSV file that
# contains the perpendicular distance from the curve and the normalized distance along the curve, for each annotated cell.
# this output is plotted in R


WHITE = (255, 255, 255)
RED = (0, 0, 255)
CELL_RADIUS = 3

draw_line_data = []
click_cells_data = []

def draw_line_mouse_callback(image, window_name, event, x, y, flags, param):
    if event == cv2.EVENT_LBUTTONDOWN:
        if len(draw_line_data) > 0:
            prev_pt = draw_line_data[-1]
            image = cv2.line(image, prev_pt, (x, y), WHITE, 3)
            cv2.imshow(window_name, image)
        draw_line_data.append((x, y))

def line(path, line_path):
    image = cv2.imread(path)
    cv2.namedWindow("image", cv2.WINDOW_AUTOSIZE)
    cv2.setMouseCallback("image", partial(draw_line_mouse_callback, image, "image"))
    cv2.imshow('image', image)
    cv2.waitKey(0)
    with open(line_path, 'w') as f:
        json.dump(draw_line_data, f)
    cv2.destroyAllWindows()

def click_cells_mouse_callback(image, window_name, event, x, y, flags, param):
    if event == cv2.EVENT_LBUTTONDOWN:
        image = cv2.circle(image, (x, y), CELL_RADIUS, RED, 3)
        cv2.imshow(window_name, image)
        click_cells_data.append((x, y))

def cells(path, cells_path):
    image = cv2.imread(path)
    cv2.namedWindow("image", cv2.WINDOW_AUTOSIZE)
    cv2.setMouseCallback("image", partial(click_cells_mouse_callback, image, "image"))
    cv2.imshow('image', image)
    cv2.waitKey(0)
    with open(cells_path, 'w') as f:
        json.dump(click_cells_data, f)
    cv2.destroyAllWindows()

def project_point_onto_line(point, line):
    """Returns:
        (i) the closest point on the line segment
        (ii) the length along the line segment to the closest point
    """
    # point = (x, y)
    # line = ((x1, y1), (x2, y2))
    c = np.array(point)
    u = np.array(line[0])
    v = np.array(line[1])
    a = c - u
    b = v - u
    scale = np.dot(a, b) / np.dot(b, b)
    a1 = scale * b
    a2 = a - a1
    if scale > 1.0:
        # beyond the end; the closest point is v
        length_along_segment = np.linalg.norm(b)
        point_on_line = v
    elif scale < 0.0:
        # before the beginning
        length_along_segment = 0.0
        point_on_line = u
    else:
        length_along_segment = np.linalg.norm(a1)
        point_on_line = u + a1
    return (point_on_line, length_along_segment)

#assert project_point_onto_line((0, 0), ((-1, -1), (-1, 1))) == (1, 1)
#assert project_point_onto_line((0, -2), ((-1, -1), (-1, 1))) == (np.sqrt(2), 0.0)
#assert project_point_onto_line((0, 2), ((-1, -1), (-1, 1))) == (np.sqrt(2), 2.0)


def process_point(lines, point):
    """Returns:
            (i) shortest distance from this point to the 'curve'
            (ii) the normalized position of that closest point along the curve (between 0 and 1)
    """
    min_distance_so_far = float('inf')
    closest_point_so_far = None

    # the sum of the lengths of the full line segments that we have visited already
    length_along_curve = 0.0

    # the length along curve of the point that gave us the min_distance_so_far
    best_length_along_curve = None

    for line in lines:
        # line is ((x1, y1), (x2, y2))
        point_on_line, length_along_segment = project_point_onto_line(point, line)
        min_distance = np.linalg.norm(point_on_line - point)
        if min_distance < min_distance_so_far:
            min_distance_so_far = min_distance
            best_length_along_curve = length_along_curve + length_along_segment
            closest_point_so_far = point_on_line
        length_along_curve += np.linalg.norm(np.array(line[1]) - np.array(line[0]))

    normalized_position = best_length_along_curve / length_along_curve
    assert normalized_position >= 0.0
    assert normalized_position <= 1.0
    return closest_point_so_far, normalized_position



def analyze(image_path, line_path, cells_path, output_csv_path, output_image_path, plot_path):
    print("analyzing...")

    with open(line_path, 'r') as f:
        line_points = json.load(f)

    if cells_path.endswith('.json'):
        with open(cells_path, 'r') as f:
            cells = json.load(f)
    elif cells_path.endswith('.tsv'):
        df = pd.read_csv(cells_path, sep='\t', names=['Type', 'Slice', 'X', 'Y', 'Value', 'C-pos', 'Z-pos', 'T-pos', 'X(pixel)', 'Y(pixel)', 'Z(pixel)'])
        x = df['X']
        y = df['Y']
        cells = list(zip(x, y))
        print(cells)
    print(line_points)

    # list of tuples
    lines = []
    prev_point = None
    for point in line_points:
        if prev_point is not None:
            lines.append((prev_point, point))
        prev_point = point

    print(f"lines: {lines}")
    print(f"cells: {cells}")

    image = cv2.imread(image_path)
    cv2.namedWindow("image", cv2.WINDOW_AUTOSIZE)
    cv2.imshow('image', image)

    # draw the curve
    for line in lines:
        (pt1, pt2) = line
        image = cv2.line(image, pt1, pt2, WHITE, 3)
    
    # draw closest point lines and populate data frame
    xs = []
    ys = []
    distances = []
    positions = []
    for cell in cells:
        closest_point, position = process_point(lines, cell)
        distance = np.linalg.norm(np.array(cell) - closest_point)
        xs.append(cell[0])
        ys.append(cell[1])
        distances.append(distance)
        positions.append(position)
        pt1 = (int(cell[0]), int(cell[1]))
        pt2 = (int(closest_point[0]), int(closest_point[1]))
        image = cv2.line(image, pt1, pt2, RED, 2)
    data = {
        "x": xs,
        "y": ys,
        "distance": distances,
        "normalized_position": positions
    }
    df = pd.DataFrame(data)
    print(df)
    df.to_csv(output_csv_path, sep="\t")

    cv2.imshow("image", image)
    cv2.waitKey(0)
    cv2.imwrite(output_image_path, image)

    # plot
    plt.figure(figsize=(8, 2.5))
    plt.scatter(positions, distances, color="green", s=10)
    plt.xlabel("normalized zonal position")
    plt.ylabel("distance (pixels)")
    plt.gca().set_xlim(0, 1)
    plt.tight_layout()
    plt.savefig(plot_path)

if len(sys.argv) < 2:
    print("Usage: python show.py <line|cells|analyze> ...") 
    exit(1)

if sys.argv[1] == "line":
    if len(sys.argv) != 4:
        print("Usage: python show.py line <input_image.png> <line.json>") 
        exit(1)
    (_, _, image_path, output_line_path) = sys.argv
    line(image_path, output_line_path)
elif sys.argv[1] == "cells":
    if len(sys.argv) != 4:
        print("Usage: python show.py cells <input_image.png> <cells.json>") 
        exit(1)
    (_, _, image_path, output_cells_path) = sys.argv
    cells(image_path, output_cells_path)
elif sys.argv[1] == "analyze":
    if len(sys.argv) != 8:
        print("Usage: python show.py analyze <input_image.png> <line.json> <cells.json> <output.tsv> <output_image.png> <plot.png>") 
        exit(1)
    (_, _, image_path, line_path, cells_path, output_csv_path, output_image_path, plot_path) = sys.argv
    analyze(image_path, line_path, cells_path, output_csv_path, output_image_path, plot_path)
else:
    print("Usage: python show.py <cells|line|analyze> ...") 
    exit(1)


