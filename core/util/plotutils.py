import numpy as np
import cv2
import os


def make_movie(image_folder, video_name):
    images = sorted([img for img in os.listdir(image_folder) if img.endswith(".png")])
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, _ = frame.shape

    video = cv2.VideoWriter(video_name, 0, 10, (width, height))

    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))

    cv2.destroyAllWindows()
    video.release()


def is_array(var):
    return hasattr(var, "__iter__") and not isinstance(var, str)


def create_dat(data, headers, fname, types=None):
    # Convert the data from list to numpy array
    if isinstance(data, list):
        data = np.array(data).T

    # Convert the data from 1D array to 2D array
    if isinstance(data, np.ndarray):
        if data.ndim == 1:
            data = data.reshape(-1, 1)

    # Handle headers where a string is passed instead of an array
    if not is_array(headers):
        if "{" in headers and "}" in headers:
            assert isinstance(data, np.ndarray)
            headers = np.repeat(headers, data.shape[1])
        else:
            headers = [headers]

    # Handle types where a class is passed instead of an array
    if not types is None:
        if not is_array(types):
            types = [types]
        assert len(headers) == len(
            types
        ), "the number of types must equal the number of headers"

    output = []

    for i, header in enumerate(headers):
        if isinstance(data, dict):
            column = data[header]
        elif isinstance(data, np.ndarray):
            headers[i] = header.format(i)
            column = data[:, i]
        else:
            raise TypeError("data must be given as a dictionary, list or numpy array")

        if types is None:
            column = np.array(column)
        else:
            column = np.array(column, dtype=types[i])

        output.append(column)

    output = np.array(output, dtype=str).T

    if fname[-4:] != ".dat":
        fname += ".dat"

    os.makedirs(os.path.dirname(fname), exist_ok=True)
    with open(fname, "w") as fmesh:
        fmesh.write(" ".join(headers) + "\n")

        for line in output:
            fmesh.write(" ".join(line) + "\n")
