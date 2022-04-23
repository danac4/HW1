import math
import sys


def euclidean_norm(point):
    x = 0
    for i in range(len(point)):
        x += math.pow(point[i], 2)
    return math.sqrt(x)


class Cluster:
    def __init__(self, centroid, new_points, cnt):
        self.centroid = centroid
        self.cnt = cnt
        self.new_points = new_points

    def format(self):
        # formatted = (["{:.4f}".format(i) for i in self.centroid])
        # return ",".join(formatted)
        tmp = [round(coord, 4) for coord in self.centroid]
        formatted = ["{:.4f}".format(i) for i in tmp]
        return ",".join(formatted)

    def append(self, point):
        for i in range(len(point)):
            self.new_points[i] += point[i]
        self.cnt += 1

    def centroid_update(self):
        if self.cnt == 0:
            return True
        new_centroid = [self.new_points[i] / self.cnt for i in range(len(self.centroid))]
        epsilon_diff = self.epsilon_check(new_centroid)
        self.centroid = new_centroid
        self.new_points = [0 for i in range(len(self.centroid))]
        self.cnt = 0
        return epsilon_diff

    def epsilon_check(self, new_centroid):
        diff = (euclidean_norm(self.centroid) - euclidean_norm(new_centroid))
        if abs(diff) >= 0.001:
            return False
        return True


def clusters_initialize(k, dim, input_points):
    points_mtx = [0 for i in range(dim)]
    clusters = []
    for i in range(k):
        pmCopy = list(points_mtx)
        clusters.append(Cluster(input_points[i], pmCopy, 0))
    return clusters


def distance(point1, point2):
    dis = 0
    for i in range(len(point1)):
        x = (point1[i] - point2[i])
        dis += math.pow(x, 2)
    res = math.sqrt(dis)
    return res


def min_distance(point, clusters):
    min_val = distance(point, clusters[0].centroid)
    min_index = 0
    for i in range(1, len(clusters)):
        curr = distance(point, clusters[i].centroid)
        if curr < min_val:
            min_index = i
            min_val = curr
    return min_index


def main_validity_check():
    assert 3 < len(sys.argv) < 6, 'Invalid Input!'
    k = validity_check_k()
    max_iter = validity_check_iter()
    points, dim = validity_check_file()
    assert k <= len(points), 'Invalid Input!'
    return k, max_iter, points, dim


def validity_check_iter():
    max_iter = 200
    if len(sys.argv) > 4:
        assert int(sys.argv[2]), 'Invalid Input!'
        max_iter = int(sys.argv[2])
        assert max_iter > 0, 'Invalid Input!'
    return max_iter


def validity_check_file():
    input_points = []
    file = sys.argv[3]
    if len(sys.argv) < 5:
        file = sys.argv[2]
    try:
        file_process = open(file, "r")
        for line in file_process:
            vector1 = line.strip()
            assert len(line) > 0, 'Invalid Input!'
            vector1 = vector1.split(",")
            try:
                vector2 = [float(num) for num in vector1]
            except ValueError:
                print('Invalid Input!')
                file_process.close()
                sys.exit(1)
            input_points.append(vector2)
        file_process.close()
        dim = len(input_points[0])
    except IOError:
        print('Invalid Input!')
        sys.exit(1)
    return input_points, dim


def validity_check_k():
    assert int(sys.argv[1]), 'Invalid Input!'
    k = int(sys.argv[1])
    assert k > 1, 'Invalid Input!'
    return k


def clusters_update(clusters):
    res = True
    for cluster in clusters:
        epsilon_indicator = cluster.centroid_update()
        res = res and epsilon_indicator
    return res


def write_to_output(result):
    output = sys.argv[-1]
    try:
        file = open(output, "a")
        for vector in result:
            file.write(vector + "\n")
        file.close()
    except IOError:
        print('Invalid Input!')
        sys.exit(1)


def kmeans():
    k, max_iter, points, dim = main_validity_check()
    clusters = clusters_initialize(k, dim, points)
    i, epsilon_indicator = 0, False
    while i < max_iter and not epsilon_indicator:
        for point in points:
            index = min_distance(point, clusters)
            clusters[index].append(point)
        epsilon_indicator = clusters_update(clusters)
        i += 1
    result = [cluster.format() for cluster in clusters]
    write_to_output(result)


if __name__ == '__main__':
    kmeans()




