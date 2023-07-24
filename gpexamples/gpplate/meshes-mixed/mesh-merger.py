import numpy as np

for topmesh in ['0', '1', '2', '3']:
    for botmesh in ['0', '1', '2', '3']:

        with open('plate_r' + topmesh + '.msh', 'r') as msh:
            lines = msh.readlines()

            # Extract the Nodes and Elements blocks from the gmsh file
            nlines = lines[lines.index('$Nodes\n')+2:lines.index('$EndNodes\n')]
            elines = lines[lines.index('$Elements\n')+2:lines.index('$EndElements\n')]

            # Split the node info
            node_ids = np.genfromtxt(nlines, dtype=int, ndmin=2)[:,0]
            coords = np.genfromtxt(nlines, dtype=float, ndmin=2)[:,1:]

            top_nodes = node_ids[np.logical_and(10000 < node_ids, node_ids < 20000)]
            top_coords = coords[np.logical_and(10000 < node_ids, node_ids < 20000)]
            tmid_nodes = node_ids[np.logical_and(20000 < node_ids, node_ids < 30000)]
            tmid_coords = coords[np.logical_and(20000 < node_ids, node_ids < 30000)]

            # Split the element info
            elem_ids = np.genfromtxt(elines, dtype=int, ndmin=2)[:,0]
            elem_info = np.genfromtxt(elines, dtype=int, ndmin=2)[:,1:5]
            inodes = np.genfromtxt(elines, dtype=int, ndmin=2)[:,5:]

            top_elems = elem_ids[np.logical_and(np.any(10000 < inodes, axis=1), np.any(inodes < 20000, axis=1))]
            top_info = elem_info[np.logical_and(np.any(10000 < inodes, axis=1), np.any(inodes < 20000, axis=1))]
            top_inodes = inodes[np.logical_and(np.any(10000 < inodes, axis=1), np.any(inodes < 20000, axis=1))]

        with open('plate_r' + botmesh + '.msh', 'r') as msh:
            lines = msh.readlines()

            # Extract the Nodes and Elements blocks from the gmsh file
            nlines = lines[lines.index('$Nodes\n')+2:lines.index('$EndNodes\n')]
            elines = lines[lines.index('$Elements\n')+2:lines.index('$EndElements\n')]

            # Split the node info
            node_ids = np.genfromtxt(nlines, dtype=int, ndmin=2)[:,0]
            coords = np.genfromtxt(nlines, dtype=float, ndmin=2)[:,1:]

            bot_nodes = node_ids[np.logical_and(30000 < node_ids, node_ids < 40000)]
            bot_coords = coords[np.logical_and(30000 < node_ids, node_ids < 40000)]
            bmid_nodes = node_ids[np.logical_and(20000 < node_ids, node_ids < 30000)]
            bmid_coords = coords[np.logical_and(20000 < node_ids, node_ids < 30000)]

            # Split the element info
            elem_ids = np.genfromtxt(elines, dtype=int, ndmin=2)[:,0]
            elem_info = np.genfromtxt(elines, dtype=int, ndmin=2)[:,1:5]
            inodes = np.genfromtxt(elines, dtype=int, ndmin=2)[:,5:]

            bot_elems = elem_ids[np.logical_and(np.any(30000 < inodes, axis=1), np.any(inodes < 40000, axis=1))]
            bot_info = elem_info[np.logical_and(np.any(30000 < inodes, axis=1), np.any(inodes < 40000, axis=1))]
            bot_inodes = inodes[np.logical_and(np.any(30000 < inodes, axis=1), np.any(inodes < 40000, axis=1))]

        if len(tmid_nodes) < len(bmid_nodes):
            cfmapping = {}
            for cnode, ccoords in zip(tmid_nodes, tmid_coords):
                for fnode, fcoords in zip(bmid_nodes, bmid_coords):
                    if np.isclose(ccoords[0], fcoords[0]):
                        cfmapping[cnode] = fnode

            if len(cfmapping) != len(tmid_nodes):
                raise ValueError('not all coarse nodes could be mapped to a fine node')

            for i in range(top_inodes.shape[0]):
                for j in range(top_inodes.shape[1]):
                    if top_inodes[i,j] in cfmapping:
                        top_inodes[i,j] = cfmapping[top_inodes[i,j]]

            mid_nodes = bmid_nodes
            mid_coords = bmid_coords

        else:
            cfmapping = {}
            for cnode, ccoords in zip(bmid_nodes, bmid_coords):
                for fnode, fcoords in zip(tmid_nodes, tmid_coords):
                    if np.isclose(ccoords[0], fcoords[0]):
                        cfmapping[cnode] = fnode

            if len(cfmapping) != len(bmid_nodes):
                raise ValueError('not all coarse nodes could be mapped to a fine node')

            for i in range(bot_inodes.shape[0]):
                for j in range(bot_inodes.shape[1]):
                    if bot_inodes[i,j] in cfmapping:
                        bot_inodes[i,j] = cfmapping[bot_inodes[i,j]]

            mid_nodes = tmid_nodes
            mid_coords = tmid_coords

        with open('plate_r' + topmesh + botmesh + '.msh', 'w') as msh:
            msh.write('$MeshFormat\n')
            msh.write('2.2 0 8\n')
            msh.write('$EndMeshFormat\n')
            msh.write('$Nodes\n')
            msh.write(str(len(top_nodes)+len(mid_nodes)+len(bot_nodes))+'\n')

            for node_id, coord in zip(top_nodes, top_coords):
                msh.write(str(node_id) + ' ' + ' '.join(map(str, coord)) + '\n')
            msh.write('\n')

            for node_id, coord in zip(mid_nodes, mid_coords):
                msh.write(str(node_id) + ' ' + ' '.join(map(str, coord)) + '\n')
            msh.write('\n')

            for node_id, coord in zip(bot_nodes, bot_coords):
                msh.write(str(node_id) + ' ' + ' '.join(map(str, coord)) + '\n')

            msh.write('$EndNodes\n')
            msh.write('$Elements\n')
            msh.write(str(len(top_elems)+len(bot_elems))+'\n')

            for elem_id, elem_inf, inode in zip(top_elems, top_info, top_inodes):
                msh.write(str(elem_id) + ' ' + ' '.join(map(str, elem_inf)) + ' ' + ' '.join(map(str, inode)) + '\n')
            msh.write('\n')

            for elem_id, elem_inf, inode in zip(bot_elems, bot_info, bot_inodes):
                msh.write(str(elem_id) + ' ' + ' '.join(map(str, elem_inf)) + ' ' + ' '.join(map(str, inode)) + '\n')

            msh.write('$EndElements\n')
