import networkx as nx
import graphviz as gv
import matplotlib.pyplot as plt

if __name__ == '__main__':
    print("Enter file with all EMU reactions: (default: all_emu)")
    all_emu_file = input()
    if not all_emu_file:
        all_emu_file = 'all_emu'

    print("Enter file with input EMUs: (default: input_emu)")
    input_emu_file = input()
    if not input_emu_file:
        input_emu_file = 'input_emu'
    input_emu = []
    with open(input_emu_file) as input_handler:
        for line in input_handler:
            input_emu.append(line[0:-1])

    print("Enter file with measured EMUs: (default: measured_emu)")
    measured_emu_file = input()
    if not measured_emu_file:
        measured_emu_file = 'measured_emu'

    measured_emu = []
    with open(measured_emu_file) as meas_handler:
        for line in meas_handler:
            measured_emu.append(line[0:-1])

    G = nx.DiGraph()
    with open(all_emu_file) as all_emu:
        for line in all_emu:
            edge = line.split(' = ')
            print(line)
            print(edge[0])
            print(edge[1][0:-1])
            print()
            if '+' in edge[0]:
                substrates = edge[0].split(' + ')
                G.add_edge(substrates[0], edge[0], color='red')
                G.add_edge(substrates[1], edge[0], color='red')
                input_emu.append(edge[0])

            G.add_edge(edge[0], edge[1][0:-1], color='black')

    colors = []
    for node in G:
        if node in measured_emu:
            colors.append('red')
        elif node in input_emu:
            colors.append('orange')
        elif G.in_degree()[node] == 0:
            colors.append('orange')
        else:
            colors.append('blue')

    edge_colors = [G[u][v]['color'] for u,v in G.edges()]

    nx.draw(G, node_color=colors, edge_color=edge_colors, with_labels=True, pos=nx.drawing.nx_agraph.graphviz_layout(G))
    plt.show()

