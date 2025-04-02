import numpy as np
import re
from ast import literal_eval
from copy import deepcopy


def parse_file(fname):
    fileraw = open(fname, "r").read()

    filestr = uncomment_file(fileraw)

    data = {}

    i = 0
    sp = filestr.split(";")

    while i < len(sp):
        line = re.sub(r"=\s*{", "={", sp[i])

        if "={" in line:
            key = line.split("={")[0]
            key = key.strip()

            newline = "={".join(line.split("={")[1:])
            data[key], i, sp = read_level(newline, i, sp)
        elif "=" in line:
            [key, value] = line.split("=")
            key = key.strip()
            value = value.strip()

            if value.startswith("[") and value.endswith("]"):
                data[key] = parse_list(value)
            else:
                data[key] = try_literal_eval(value)
        elif line != "":
            raise RuntimeError("Unable to parse: %s" % line)

        i = i + 1

    return data


def read_level(line, i, sp):
    subdata = {}

    while True:
        if "={" in line:
            key = line.split("={")[0]
            key = key.strip()

            newline = "={".join(line.split("={")[1:])
            subdata[key], i, sp = read_level(newline, i, sp)
        elif "=" in line:
            [key, value] = line.split("=")
            key = key.strip()
            value = value.strip()

            if value.startswith("[") and value.endswith("]"):
                subdata[key] = parse_list(value)
            else:
                subdata[key] = try_literal_eval(value)
        elif "}" in line:
            return subdata, i, sp

        elif line != "":
            raise RuntimeError("Unable to parse: %s" % line)

        i = i + 1

        if i == len(sp):
            raise RuntimeError(
                "EOF reached while parsing an input block. Did you forget to close a bracket?"
            )

        line = re.sub(r"=\s*{", "={", sp[i])


def uncomment_file(fileraw):
    filestr = ""
    comment_mode = False

    # Go through all lines in the raw file
    for line in fileraw.split("\n"):
        # Check if we are in comment mode
        if comment_mode:
            # If so, try to find '*/' and remove only the part before it
            end = line.find("*/")
            if end >= 0:
                comment_mode = False
                line = line[end + len("*/") :]

            # If comment_mode is still enabled, don't include anything
            else:
                line = ""

        if not comment_mode:
            # If we are not in comment mode, remove all full comments from the line
            line = uncomment_line(line)

            # If there is a '/*' in the line, enable comment mode, and remove the part after '/*' from the line
            start = line.find("/*")
            if start >= 0:
                comment_mode = True
                line = line[:start]

        # Add the line to the file string
        filestr += line

    return filestr


def uncomment_line(line):
    clean_line = line

    # Remove all comments from the line (assuming that comment_mode is False)
    while True:
        # If the first identifier is a '//', remove everything after it
        start_oneline = clean_line.find("//")
        start_block = clean_line.find("/*")
        if start_oneline >= 0:
            if start_oneline < start_block or start_block < 0:
                clean_line = clean_line[:start_oneline]

        # Remove everything in the first one-line block comment that is found
        start = clean_line.find("/*")
        end = clean_line.find("*/")
        if 0 <= start < end and end >= 0:
            clean_line = clean_line[:start] + clean_line[end + len("*/") :]
        else:
            # Exit if no comments are left
            break

    return clean_line


def parse_list(value):
    if not value.startswith("[") or not value.endswith("]"):
        raise ValueError("Not a valid list")

    lst = []
    val = ""
    depth = 0
    for s in value:
        if s == "[":
            depth += 1

        if s not in [",", "[", "]"] or depth > 1:
            val += s
        elif depth == 1 and s in [",", "]"]:
            val = val.strip()
            if val.startswith("[") and val.endswith("]"):
                val = parse_list(val)
            else:
                val = try_literal_eval(val)
            lst.append(val)
            val = ""

        if s == "]":
            depth -= 1

    return lst


def try_literal_eval(val):
    try:
        return literal_eval(val)
    except:
        return val


def props_to_file(props, fname):
    string = props_to_string(props)

    file = open(fname, "w")
    file.write(string)
    file.close()


def indent_to_string(depth):
    return " " * (2 * depth)


def props_to_string(dic, depth=0):
    string = ""

    if depth > 0:
        string += indent_to_string(depth - 1)
        string += "{\n"

    *_, lastkey = dic.keys()

    for key, value in dic.items():
        string += indent_to_string(depth)

        if isinstance(value, dict):
            string += str(key) + " =\n"

            depth += 1
            string += props_to_string(value, depth)
            depth -= 1

            if key != lastkey:
                string += "\n"

        elif isinstance(value, list):
            string += str(key) + " = "
            string += list_to_string(value)

        elif isinstance(value, str):
            string += str(key) + ' = "' + value + '";\n'

        elif isinstance(value, bool):
            string += str(key) + " = " + str(value).lower() + ";\n"

        else:
            string += str(key) + " = " + str(value) + ";\n"

    if depth > 0:
        string += indent_to_string(depth - 1)
        string += "};\n"

    return string


def list_to_string(lst):
    string = "[ "

    for i, value in enumerate(lst):
        if isinstance(value, str):
            string += '"' + value + '"'
        else:
            string += str(value)

        if i < len(lst) - 1:
            string += ", "
        else:
            string += " ];\n"

    return string


def soft_cast(value, typ):
    # This function attempts to convert value to typ
    # If this conversion fails, it returns the original value
    try:
        return typ(value)
    except:
        return value


def evaluate(value, coords, rank, extra_dict=None):
    # This function does a string evaluation of value, if possible
    if type(value) is str:
        eval_dict = get_eval_dict(coords, rank, extra_dict)
        return eval(value, {}, eval_dict)
    else:
        return value


def get_eval_dict(coords, rank, extra_dict=None):
    # This function builds a dictionary with the x,y,z coordinates of coords
    eval_dict = {"x": coords[0]}

    if rank >= 2:
        eval_dict.update({"y": coords[1]})
    if rank == 3:
        eval_dict.update({"z": coords[2]})

    # Add the extra dictionary if applicable
    if extra_dict is not None:
        eval_dict.update(extra_dict)

    return eval_dict


def get_core_eval_dict():
    core_eval_dict = {
        "exp": np.exp,
        "sin": np.sin,
        "cos": np.cos,
        "tan": np.tan,
        "pi": np.pi,
        "sqrt": np.sqrt,
    }

    return core_eval_dict


def check_dict(obj, dic, keys=[]):
    name = obj.__class__.__name__
    if not isinstance(dic, dict):
        raise ValueError("Argument in {} must be a dict".format(name))
    for key in keys:
        if key not in dic.keys():
            raise ValueError("Argument in {} must contain '{}' key".format(name, key))


def check_list(obj, lst):
    name = obj.__class__.__name__
    if not isinstance(lst, list):
        raise ValueError("Argument in {} must be a list".format(name))


def check_value(obj, val, options=[]):
    name = obj.__class__.__name__
    if val is None:
        raise ValueError("Argument in {} cannot be 'None'".format(name))
    if len(options) > 0 and val not in options:
        raise ValueError("Argument in {} must be one of {}".format(name, options))


def split_off_type(props):
    cprops = deepcopy(props)
    typ = cprops.pop("type")
    return typ, cprops


def get_recursive(dic, keys):
    sub = dic
    for key in keys:
        sub = sub[key]
    return sub


def set_recursive(dic, keys, value):
    sub = dic
    for i, key in enumerate(keys):
        if i == len(keys) - 1:
            sub[key] = value
        else:
            if key not in sub:
                sub[key] = {}
            sub = sub[key]


def get_attr_recursive(obj, keys):
    attr = obj
    for key in keys:
        attr = getattr(attr, key)
    return attr


def set_attr_recursive(obj, keys, value):
    attr = obj
    for key in keys[:-1]:
        attr = getattr(attr, key)
    setattr(attr, keys[-1], value)


def call_attr_recursive(obj, keys, value):
    attr = obj
    for key in keys:
        attr = getattr(attr, key)
    attr(value)


def set_or_call_attr_recursive(obj, keys, value):
    attr = obj
    for key in keys[:-1]:
        attr = getattr(attr, key)
    if callable(getattr(attr, keys[-1])):
        attr = getattr(attr, keys[-1])
        attr(value)
    else:
        setattr(attr, keys[-1], value)


def split_key(key):
    keys = key.split(".")
    for i, k in enumerate(keys):
        if k.lstrip("+-").isdigit():
            keys[i] = int(k)
    return keys
