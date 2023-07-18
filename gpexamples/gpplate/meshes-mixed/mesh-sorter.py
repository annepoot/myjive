import numpy as np

for fname in ['plate_r0.msh', 'plate_r1.msh', 'plate_r2.msh']:
    with open(fname, 'r') as msh:
        lines = msh.readlines()

        # Extract the Nodes and Elements blocks from the gmsh file
        nlines = lines[lines.index('$Nodes\n')+2:lines.index('$EndNodes\n')]
        elines = lines[lines.index('$Elements\n')+2:lines.index('$EndElements\n')]

        # Split the node info
        node_ids = np.genfromtxt(nlines, dtype=int, ndmin=2)[:,0]
        coords = np.genfromtxt(nlines, dtype=float, ndmin=2)[:,1:]

        top_nodes = []
        mid_nodes = []
        bot_nodes = []
        top_coords = []
        mid_coords = []
        bot_coords = []

        for i, (node_id, coord) in enumerate(zip(node_ids, coords[:,1])):
            if np.isclose(coord, 1):
                mid_nodes.append(node_id)
                mid_coords.append(coords[i,:])
            elif coord > 1:
                top_nodes.append(node_id)
                top_coords.append(coords[i,:])
            elif coord < 1:
                bot_nodes.append(node_id)
                bot_coords.append(coords[i,:])

        top_nodes = np.array(top_nodes)
        mid_nodes = np.array(mid_nodes)
        bot_nodes = np.array(bot_nodes)

        top_coords = np.array(top_coords)
        mid_coords = np.array(mid_coords)
        bot_coords = np.array(bot_coords)

        # Split the element info
        elem_ids = np.genfromtxt(elines, dtype=int, ndmin=2)[:,0]
        elem_info = np.genfromtxt(elines, dtype=int, ndmin=2)[:,1:5]
        inodes = np.genfromtxt(elines, dtype=int, ndmin=2)[:,5:]

        if np.all(top_nodes < 1000) and np.all(mid_nodes < 1000) and np.all(bot_nodes < 1000):
            reindex = True
        elif np.all(1000 < top_nodes) and np.all(top_nodes < 2000) and np.all(2000 < mid_nodes) and np.all(mid_nodes < 3000) and np.all(3000 < bot_nodes) and np.all(bot_nodes < 4000):
            reindex = False
        else:
            raise ValueError('Cannot decide if node ids should be updated or not.')

        if reindex:
            for i in range(inodes.shape[0]):
                for j in range(inodes.shape[1]):
                    if inodes[i,j] in top_nodes:
                        inodes[i,j] += 1000
                    elif inodes[i,j] in mid_nodes:
                        inodes[i,j] += 2000
                    elif inodes[i,j] in bot_nodes:
                        inodes[i,j] += 3000
                    else:
                        raise ValueError('node was not found in top, mid or bottom nodes')

            top_nodes += 1000
            mid_nodes += 2000
            bot_nodes += 3000

            elem_ids[np.logical_and(np.any(1000 < inodes, axis=1), np.any(inodes < 2000, axis=1))] += 1000
            elem_ids[np.logical_and(np.any(3000 < inodes, axis=1), np.any(inodes < 4000, axis=1))] += 3000

    with open(fname, 'w') as msh:
        for line in lines[:lines.index('$Nodes\n')+2]:
            msh.write(line)

        for node_id, coord in zip(top_nodes, top_coords):
            msh.write(str(node_id) + ' ' + ' '.join(map(str, coord)) + '\n')
        msh.write('\n')

        for node_id, coord in zip(mid_nodes, mid_coords):
            msh.write(str(node_id) + ' ' + ' '.join(map(str, coord)) + '\n')
        msh.write('\n')

        for node_id, coord in zip(bot_nodes, bot_coords):
            msh.write(str(node_id) + ' ' + ' '.join(map(str, coord)) + '\n')
        msh.write('\n')

        for line in lines[lines.index('$EndNodes\n'):lines.index('$Elements\n')+2]:
            msh.write(line)

        for elem_id, elem_inf, inode in zip(elem_ids, elem_info, inodes):
            msh.write(str(elem_id) + ' ' + ' '.join(map(str, elem_inf)) + ' ' + ' '.join(map(str, inode)) + '\n')

        for line in lines[lines.index('$EndElements\n'):]:
            msh.write(line)
