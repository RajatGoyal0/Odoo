// ProceduralTerrain.cpp
// Production-ready, single-file procedural terrain generator for open-world adventure games.
// C++23, depends on SDL2.
// Compile with: g++ -std=c++23 ProceduralTerrain.cpp -lSDL2 -o ProceduralTerrain

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include <queue>
#include <cassert>
#include <SDL2/SDL.h>

// Configurable constants
constexpr int MAP_WIDTH = 128;
constexpr int MAP_HEIGHT = 128;
constexpr int TILE_SIZE = 5;           // Pixels per tile in visualizer
constexpr float MAX_ELEVATION = 50.f;
constexpr int NUM_POIS = 5;            // Key Points of Interest

// Tile type enum
enum class TileType { Empty, Walkable, POI };

// Tile struct
struct Tile {
    float elevation = 0.f;
    TileType type = TileType::Empty;
};

// Noise generator (simple Perlin-like)
class Noise {
    std::vector<int> perm;
    float fade(float t) const { return t * t * t * (t * (t * 6 - 15) + 10); }
    float grad(int hash, float x, float y) const {
        int h = hash & 7; // 8 directions
        float u = h < 4 ? x : y;
        float v = h < 4 ? y : x;
        return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
    }
public:
    Noise(unsigned seed) {
        perm.resize(512);
        std::iota(perm.begin(), perm.begin() + 256, 0);
        std::shuffle(perm.begin(), perm.begin() + 256, std::mt19937(seed));
        for (int i = 0; i < 256; ++i) perm[256 + i] = perm[i];
    }
    float get(float x, float y) const {
        int X = (int)std::floor(x) & 255;
        int Y = (int)std::floor(y) & 255;
        x -= std::floor(x);
        y -= std::floor(y);
        float u = fade(x), v = fade(y);
        int A = perm[X] + Y, B = perm[X + 1] + Y;
        return std::lerp(
            std::lerp(grad(perm[A], x, y), grad(perm[B], x - 1, y), u),
            std::lerp(grad(perm[A + 1], x, y - 1), grad(perm[B + 1], x - 1, y - 1), u),
            v
        );
    }
};

// Main generator
class ProceduralLevel {
    std::vector<std::vector<Tile>> map;
    std::vector<std::pair<int,int>> pois;
    Noise noise;

    bool isValid(int x, int y) const {
        return x >= 0 && y >= 0 && x < MAP_WIDTH && y < MAP_HEIGHT;
    }

    // BFS to ensure paths between POIs
    void ensurePaths() {
        for (size_t i = 1; i < pois.size(); ++i) {
            auto [sx, sy] = pois[i - 1];
            auto [ex, ey] = pois[i];
            bfsPath(sx, sy, ex, ey);
        }
    }

    void bfsPath(int sx, int sy, int ex, int ey) {
        std::queue<std::pair<int,int>> q;
        std::vector<std::vector<bool>> visited(MAP_WIDTH, std::vector<bool>(MAP_HEIGHT, false));
        std::vector<std::vector<std::pair<int,int>>> prev(MAP_WIDTH, std::vector<std::pair<int,int>>(MAP_HEIGHT, {-1,-1}));
        q.emplace(sx, sy);
        visited[sx][sy] = true;
        const int dx[] = {1,-1,0,0}, dy[] = {0,0,1,-1};

        while (!q.empty()) {
            auto [x, y] = q.front(); q.pop();
            if (x == ex && y == ey) break;
            for (int d=0; d<4; ++d) {
                int nx = x+dx[d], ny=y+dy[d];
                if (isValid(nx,ny) && !visited[nx][ny]) {
                    visited[nx][ny]=true;
                    prev[nx][ny]={x,y};
                    q.emplace(nx,ny);
                }
            }
        }
        // reconstruct path
        for (int x=ex, y=ey; x!=-1 && y!=-1; std::tie(x,y)=prev[x][y]) {
            map[x][y].type = TileType::Walkable;
        }
    }

public:
    ProceduralLevel(unsigned seed=42) : noise(seed) {
        map.resize(MAP_WIDTH, std::vector<Tile>(MAP_HEIGHT));
    }

    void generate() {
        // Generate terrain using noise
        for (int x=0; x<MAP_WIDTH; ++x) {
            for (int y=0; y<MAP_HEIGHT; ++y) {
                float n = noise.get(x * 0.1f, y * 0.1f);
                n = (n+1.f)/2.f; // normalize 0..1
                map[x][y].elevation = n * MAX_ELEVATION;
                map[x][y].type = (n>0.3f) ? TileType::Walkable : TileType::Empty;
            }
        }
        // Place POIs
        std::mt19937 rng(123);
        std::uniform_int_distribution<> distX(10, MAP_WIDTH-10), distY(10, MAP_HEIGHT-10);
        pois.clear();
        for (int i=0;i<NUM_POIS;++i) {
            int x=distX(rng), y=distY(rng);
            pois.emplace_back(x,y);
            map[x][y].type=TileType::POI;
        }
        // Ensure paths
        ensurePaths();
    }

    void save(const std::string& filename) {
        std::ofstream out(filename, std::ios::binary);
        for (const auto& row: map)
            for (const auto& tile: row)
                out.write(reinterpret_cast<const char*>(&tile), sizeof(Tile));
    }

    void load(const std::string& filename) {
        std::ifstream in(filename, std::ios::binary);
        for (auto& row: map)
            for (auto& tile: row)
                in.read(reinterpret_cast<char*>(&tile), sizeof(Tile));
    }

    const std::vector<std::vector<Tile>>& getMap() const { return map; }
};

// Minimal SDL2 visualizer
void visualize(const ProceduralLevel& level) {
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* win = SDL_CreateWindow("Terrain Preview", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, MAP_WIDTH*TILE_SIZE, MAP_HEIGHT*TILE_SIZE, 0);
    SDL_Renderer* ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
    auto& map = level.getMap();

    bool running=true;
    while (running) {
        SDL_Event e;
        while(SDL_PollEvent(&e)) {
            if (e.type==SDL_QUIT) running=false;
        }
        SDL_SetRenderDrawColor(0,0,0,255); SDL_RenderClear(ren);
        for (int x=0;x<MAP_WIDTH;++x) {
            for (int y=0;y<MAP_HEIGHT;++y) {
                auto& t=map[x][y];
                if (t.type==TileType::POI)
                    SDL_SetRenderDrawColor(255,0,0,255);
                else if (t.type==TileType::Walkable)
                    SDL_SetRenderDrawColor(0,200,0,255);
                else
                    SDL_SetRenderDrawColor(50,50,50,255);
                SDL_Rect r{x*TILE_SIZE, y*TILE_SIZE, TILE_SIZE, TILE_SIZE};
                SDL_RenderFillRect(ren,&r);
            }
        }
        SDL_RenderPresent(ren);
        SDL_Delay(16); // ~60 FPS
    }
    SDL_DestroyRenderer(ren); SDL_DestroyWindow(win); SDL_Quit();
}

int main() {
    ProceduralLevel level;
    level.generate();
    level.save("level.dat");
    level.load("level.dat");
    visualize(level);
    std::cout << "Level generated, saved, loaded and visualized successfully.\n";
    return 0;
}
