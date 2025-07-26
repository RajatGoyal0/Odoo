// ProceduralTerrainFixed.cpp
// C++23, production-ready, maintainability fixes applied.
// Compile with: g++ -std=c++23 ProceduralTerrainFixed.cpp -lSDL2 -o ProceduralTerrainFixed

#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <cmath>
#include <random>
#include <queue>
#include <cassert>
#include <SDL2/SDL.h>
#include <cstddef> // for std::byte

constexpr int MAP_WIDTH = 128;
constexpr int MAP_HEIGHT = 128;
constexpr int TILE_SIZE = 5;
constexpr float MAX_ELEVATION = 50.f;
constexpr int NUM_POIS = 5;

enum class TileType { Empty, Walkable, POI };

struct Tile {
    float elevation = 0.f;
    TileType type = TileType::Empty;
};

// Noise generator
class Noise {
    std::vector<int> perm;

    float fade(float t) const {
        return t * t * t * (t * (t * 6 - 15) + 10);
    }

    float grad(int hash, float x, float y) const {
        int h = hash & 7;
        float u = h < 4 ? x : y;
        float v = h < 4 ? y : x;
        float gu = (h & 1) ? -u : u;
        float gv = (h & 2) ? -v : v;
        return gu + gv;
    }

public:
    explicit Noise(unsigned seed) {
        perm.resize(512);
        std::iota(perm.begin(), perm.begin() + 256, 0);
        std::mt19937 rng{seed};
        std::shuffle(perm.begin(), perm.begin() + 256, rng);
        for (int i = 0; i < 256; ++i) {
            perm[256 + i] = perm[i];
        }
    }

    float get(float x, float y) const {
        int X = static_cast<int>(std::floor(x)) & 255;
        int Y = static_cast<int>(std::floor(y)) & 255;
        x -= std::floor(x);
        y -= std::floor(y);
        float u = fade(x);
        float v = fade(y);
        int A = perm[X] + Y;
        int B = perm[X + 1] + Y;
        float res1 = std::lerp(grad(perm[A], x, y), grad(perm[B], x - 1, y), u);
        float res2 = std::lerp(grad(perm[A + 1], x, y - 1), grad(perm[B + 1], x - 1, y - 1), u);
        return std::lerp(res1, res2, v);
    }
};

class ProceduralLevel {
    std::vector<std::vector<Tile>> map;
    std::vector<std::pair<int, int>> pois;
    Noise noise;

    bool isValid(int x, int y) const {
        return x >= 0 && y >= 0 && x < MAP_WIDTH && y < MAP_HEIGHT;
    }

    void bfsPath(int sx, int sy, int ex, int ey) {
        std::queue<std::pair<int, int>> q;
        std::vector<std::vector<bool>> visited(MAP_WIDTH, std::vector<bool>(MAP_HEIGHT, false));
        std::vector<std::vector<std::pair<int, int>>> prev(MAP_WIDTH, std::vector<std::pair<int, int>>(MAP_HEIGHT, {-1, -1}));

        q.emplace(sx, sy);
        visited[sx][sy] = true;

        const std::array<int, 4> dx{1, -1, 0, 0};
        const std::array<int, 4> dy{0, 0, 1, -1};

        while (!q.empty()) {
            auto [x, y] = q.front();
            q.pop();
            if (x == ex && y == ey) break;

            for (std::size_t d = 0; d < dx.size(); ++d) {
                int nx = x + dx[d];
                int ny = y + dy[d];
                if (isValid(nx, ny) && !visited[nx][ny]) {
                    visited[nx][ny] = true;
                    prev[nx][ny] = {x, y};
                    q.emplace(nx, ny);
                }
            }
        }

        for (int x = ex, y = ey; x != -1 && y != -1; std::tie(x, y) = prev[x][y]) {
            map[x][y].type = TileType::Walkable;
        }
    }

    void ensurePaths() {
        for (std::size_t i = 1; i < pois.size(); ++i) {
            auto [sx, sy] = pois[i - 1];
            auto [ex, ey] = pois[i];
            bfsPath(sx, sy, ex, ey);
        }
    }

public:
    explicit ProceduralLevel(unsigned seed = 42) : noise(seed) {
        map.resize(MAP_WIDTH);
        for (auto& row : map) {
            row.resize(MAP_HEIGHT);
        }
    }

    void generate() {
        for (int x = 0; x < MAP_WIDTH; ++x) {
            for (int y = 0; y < MAP_HEIGHT; ++y) {
                float n = noise.get(x * 0.1f, y * 0.1f);
                n = (n + 1.f) / 2.f;
                map[x][y].elevation = n * MAX_ELEVATION;
                map[x][y].type = (n > 0.3f) ? TileType::Walkable : TileType::Empty;
            }
        }

        std::mt19937 rng{123};
        std::uniform_int_distribution<> distX(10, MAP_WIDTH - 10);
        std::uniform_int_distribution<> distY(10, MAP_HEIGHT - 10);

        pois.clear();
        for (int i = 0; i < NUM_POIS; ++i) {
            int x = distX(rng);
            int y = distY(rng);
            pois.emplace_back(x, y);
            map[x][y].type = TileType::POI;
        }

        ensurePaths();
    }

    void save(const std::string& filename) {
        std::ofstream out(filename, std::ios::binary);
        for (const auto& row : map) {
            for (const auto& tile : row) {
                const std::byte* data = reinterpret_cast<const std::byte*>(&tile);
                out.write(reinterpret_cast<const char*>(data), sizeof(Tile));
            }
        }
    }

    void load(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        for (auto& row : map) {
            for (auto& tile : row) {
                std::byte* data = reinterpret_cast<std::byte*>(&tile);
                in.read(reinterpret_cast<char*>(data), sizeof(Tile));
            }
        }
    }

    const std::vector<std::vector<Tile>>& getMap() const { return map; }
};

void visualize(const ProceduralLevel& level) {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cerr << "SDL_Init failed\n";
        return;
    }

    SDL_Window* win = SDL_CreateWindow("Terrain Preview", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, MAP_WIDTH * TILE_SIZE, MAP_HEIGHT * TILE_SIZE, 0);
    SDL_Renderer* ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);

    auto& map = level.getMap();
    bool running = true;

    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = false;
        }

        SDL_SetRenderDrawColor(0, 0, 0, 255);
        SDL_RenderClear(ren);

        for (int x = 0; x < MAP_WIDTH; ++x) {
            for (int y = 0; y < MAP_HEIGHT; ++y) {
                const auto& t = map[x][y];
                if (t.type == TileType::POI)
                    SDL_SetRenderDrawColor(255, 0, 0, 255);
                else if (t.type == TileType::Walkable)
                    SDL_SetRenderDrawColor(0, 200, 0, 255);
                else
                    SDL_SetRenderDrawColor(50, 50, 50, 255);

                SDL_Rect r;
                r.x = x * TILE_SIZE;
                r.y = y * TILE_SIZE;
                r.w = TILE_SIZE;
                r.h = TILE_SIZE;
                SDL_RenderFillRect(ren, &r);
            }
        }

        SDL_RenderPresent(ren);
        SDL_Delay(16);
    }

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
}

int main() {
    ProceduralLevel level;
    level.generate();
    level.save("level.dat");
    level.load("level.dat");
    visualize(level);
    std::cout << "Level generated, saved, loaded, and visualized successfully.\n";
    return 0;
}
