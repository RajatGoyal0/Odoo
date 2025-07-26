#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <random>
#include <array>
#include <string>
#include <cstdint>
#include <cstddef>

constexpr int MAP_WIDTH = 100;
constexpr int MAP_HEIGHT = 100;
constexpr int MAX_ELEVATION = 20;

// Modular tile structure
struct Tile {
    int elevation{};
    bool walkable{true};
    char type{'G'}; // G=Grass, W=Water, M=Mountain

    explicit Tile(int elevation = 0, bool walkable = true, char type = 'G')
        : elevation(elevation), walkable(walkable), type(type) {}
};

// Level data
class Level {
public:
    Level() : tiles(MAP_WIDTH, std::vector<Tile>(MAP_HEIGHT)) {}

    void generate() {
        std::mt19937 rng{123};
        auto distElevation = std::uniform_int_distribution<>(0, MAX_ELEVATION);
        auto distType = std::uniform_int_distribution<>(0, 2);

        for (auto& column : tiles) {
            for (auto& tile : column) {
                tile.elevation = distElevation(rng);
                tile.walkable = tile.elevation < (MAX_ELEVATION / 2);
                tile.type = ['G', 'W', 'M'][distType(rng)];
            }
        }

        placePointsOfInterest(rng);
        ensurePaths();
    }

    void save(const std::string& filename) const {
        std::ofstream out(filename, std::ios::binary);
        for (const auto& column : tiles) {
            for (const auto& tile : column) {
                const auto* data = reinterpret_cast<const std::byte*>(&tile);
                out.write(reinterpret_cast<const char*>(data), sizeof(Tile));
            }
        }
    }

    void load(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        for (auto& column : tiles) {
            for (auto& tile : column) {
                auto* data = reinterpret_cast<std::byte*>(&tile);
                in.read(reinterpret_cast<char*>(data), sizeof(Tile));
            }
        }
    }

    void debugVisualize() const {
        for (std::size_t y = 0; y < MAP_HEIGHT; ++y) {
            for (std::size_t x = 0; x < MAP_WIDTH; ++x) {
                std::cout << tiles[x][y].type;
            }
            std::cout << '\n';
        }
    }

private:
    std::vector<std::vector<Tile>> tiles;
    std::vector<std::pair<int, int>> pointsOfInterest;

    void placePointsOfInterest(std::mt19937& rng) {
        auto distX = std::uniform_int_distribution<>(10, MAP_WIDTH - 10);
        auto distY = std::uniform_int_distribution<>(10, MAP_HEIGHT - 10);

        for (int i = 0; i < 5; ++i) {
            pointsOfInterest.emplace_back(distX(rng), distY(rng));
        }
    }

    void ensurePaths() {
        static constexpr std::array dx{0, 1, 0, -1};
        static constexpr std::array dy{1, 0, -1, 0};

        for (std::size_t i = 1; i < pointsOfInterest.size(); ++i) {
            const auto& [sx, sy] = pointsOfInterest[i - 1];
            const auto& [ex, ey] = pointsOfInterest[i];

            auto visited = std::vector(MAP_WIDTH, std::vector(MAP_HEIGHT, false));
            auto prev = std::vector(MAP_WIDTH, std::vector(MAP_HEIGHT, std::pair{-1, -1}));
            std::queue<std::pair<int, int>> q;

            visited[sx][sy] = true;
            q.emplace(sx, sy);

            while (!q.empty()) {
                auto [x, y] = q.front();
                q.pop();
                if (x == ex && y == ey) break;

                for (std::size_t d = 0; d < dx.size(); ++d) {
                    int nx = x + dx[d];
                    int ny = y + dy[d];
                    if (!isValid(nx, ny) || visited[nx][ny]) continue;

                    visited[nx][ny] = true;
                    prev[nx][ny] = {x, y};
                    q.emplace(nx, ny);
                }
            }

            tracePath(prev, sx, sy, ex, ey);
        }
    }

    bool isValid(int x, int y) const {
        return x >= 0 && y >= 0 && x < MAP_WIDTH && y < MAP_HEIGHT && tiles[x][y].walkable;
    }

    void tracePath(const auto& prev, int sx, int sy, int ex, int ey) {
        int x = ex, y = ey;
        while (!(x == sx && y == sy)) {
            auto [px, py] = prev[x][y];
            if (px == -1 && py == -1) break;
            tiles[x][y].type = 'P'; // P = Path
            x = px;
            y = py;
        }
    }
};

int main() {
    Level level;
    level.generate();
    level.save("level.dat");

    Level loadedLevel;
    loadedLevel.load("level.dat");

    loadedLevel.debugVisualize();
    return 0;
}
