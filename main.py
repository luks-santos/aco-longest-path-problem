import pandas as pd
import numpy as np
import sys
import time

NUM_ANTS = 190
NUM_ITERATIONS = 75
INITIAL_PHEROMONE = 0.0001

RHO = 0.6
ALPHA = 0.7
BETA = 1
EXPLORATION = NUM_ANTS * 0.6


def read_graph_data(filename):
    df = pd.read_csv(filename, sep='\t')
    return df['origem'].to_numpy(), df['destino'].to_numpy(), df['custo'].to_numpy()


def aco(sources, destinations, costs, initial_vertex, end_vertex):
    start_time = time.time()

    pheromones = np.full(len(sources), INITIAL_PHEROMONE, dtype=np.float64)
    best_path = []
    best_path_cost = 0
    value_exploration = 0.1

    for it in range(NUM_ITERATIONS):
        progress = (it * 100) / NUM_ITERATIONS
        print(f'\r Iteração {it + 1}. Concluído {progress:.2f}% ', end='', flush=True)

        ant_paths, distance_paths, index_paths = [], [], []
        ant = 0
        while ant < NUM_ANTS:
            completed_path, invalid_path = False, False
            current_vertex = initial_vertex

            ant_path, ant_path_index = [], []
            ant_path_cost = 0

            while not completed_path:
                ant_path.append(current_vertex)

                vertex_indices = np.where(sources == current_vertex)[0].tolist()
                unvisited_vertices = [index for index in vertex_indices if destinations[index] not in ant_path]
                vertex_indices = unvisited_vertices

                if len(vertex_indices) == 0:
                    ant -= 1
                    completed_path = True
                    invalid_path = True
                else:
                    probabilities = (pheromones[vertex_indices] ** ALPHA) * (
                                costs[vertex_indices] ** (BETA + value_exploration))
                    probabilities /= sum(probabilities)

                    next_vertex_index = np.random.choice(vertex_indices, p=probabilities)

                    ant_path_cost += costs[next_vertex_index]
                    ant_path_index.append(next_vertex_index)

                    current_vertex = destinations[next_vertex_index]
                    if current_vertex == end_vertex:
                        completed_path = True
                        ant_path.append(current_vertex)

            if not invalid_path:
                if ant_path_cost > best_path_cost:
                    best_path_cost = ant_path_cost
                    best_path = ant_path

                ant_paths.append(ant_path)
                distance_paths.append(ant_path_cost)
                index_paths.append(ant_path_index)

            ant += 1

        pheromones *= (1 - RHO)

        for i, index_path in enumerate(index_paths):
            sum_phe = pheromones[index_path] + distance_paths[i] / 1000

            pheromones[index_path] = sum_phe

            equal_path_count = index_paths.count(index_path)

            if equal_path_count >= EXPLORATION:
                value_exploration = min(value_exploration * 2, 15)

    end_time = time.time()
    return best_path, best_path_cost, (end_time - start_time)


if __name__ == "__main__":

    base = "C"
    base_name = ""

    graph_data = []
    initial_vertex, end_vertex = 0, 0

    if base == "A":
        graph_data = read_graph_data('grafo1.csv')
        initial_vertex = 1
        end_vertex = 12
        base_name = "Graph B - 12 vertices and 25 edges"

    elif base == "B":
        graph_data = read_graph_data('grafo2.csv')
        initial_vertex = 1
        end_vertex = 20
        base_name = "Graph C - 20 vertices and 190 edges"

    elif base == "C":
        graph_data = read_graph_data('grafo3.csv')
        initial_vertex = 1
        end_vertex = 100
        base_name = "Graph D - 100 vertices and 8020 edges"

    path, best, time = aco(graph_data[0], graph_data[1], graph_data[2], initial_vertex, end_vertex)

    print(path, best, time)