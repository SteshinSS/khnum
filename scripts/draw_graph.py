import networkx as nx
import graphviz as gv
import matplotlib.pyplot as plt

if __name__ == '__main__':
    print("Enter file with all EMU reactions:")
    all_emu_file = input()

    print("Enter file with input EMUs:")
    input_emu_file = input()
    input_emu = []
    if input_emu_file:
        with open(input_emu_file) as input_handler:
            for line in input_handler:
                input_emu.append(line[0:-1])

    print("Enter file with measured EMUs:")
    measured_emu_file = input()
    measured_emu = []
    if measured_emu_file:
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
                G.add_edge(substrates[0], edge[0])
                G.add_edge(substrates[1], edge[0])

            G.add_edge(edge[0], edge[1][0:-1])

    colors = []
    for node in G:
        if node in measured_emu:
            colors.append('red')
        elif node in input_emu:
            colors.append('orange')
        else:
            colors.append('blue')

    nx.draw(G, node_color=colors, with_labels=True, pos=nx.drawing.nx_agraph.graphviz_layout(G))
    plt.show()

