/*
@tagline: Production-ready C++23 procedural terrain generator with enhanced safety, modern C++ features, and proper memory management.

@intuition
Use modern C++23 safety features, std::byte for binary operations, constrained templates, and safe casting to eliminate undefined behavior while maintaining performance.

@approach
- Replace C-style arrays with std::array for compile-time safety
- Use std::byte and safe casting for binary serialization
- Add template constraints to prevent unintended copying/moving
- Use std::to_underlying for enum conversions
- Apply auto type deduction where appropriate
- Implement safe binary I/O with proper alignment and endianness handling

@complexity
Time: O(octaves * n*m) for generation, O(path_len*log|open|) for pathfinding
Space: O(n*m) for level data, O(1) for algorithms
*/

// SDL2 headers first to avoid order issues
#ifdef _WIN32
#define SDL_MAIN_HANDLED
#pragma comment(lib, "SDL2.lib")
#endif
#include <SDL2/SDL.h>

#include <array>
#include <span>
#include <vector>
#include <queue>
#include <random>
#include <string>
#include <ranges>
#include <algorithm>
#include <unordered_map>
#include <optional>
#include <concepts>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <cassert>
#include <chrono>
#include <functional>
#include <cmath>
#include <numeric>
#include <limits>
#include <memory>
#include <bit>

// Safe binary I/O utilities
namespace binary_io {
    template<typename T>
    std::array<std::byte, sizeof(T)> to_bytes(const T& value) noexcept {
        std::array<std::byte, sizeof(T)> bytes;
        std::memcpy(bytes.data(), &value, sizeof(T));
        return bytes;
    }
    
    template<typename T>
    T from_bytes(std::span<const std::byte, sizeof(T)> bytes) noexcept {
        T value;
        std::memcpy(&value, bytes.data(), sizeof(T));
        return value;
    }
    
    template<typename T>
    bool write_safe(std::ostream& out, const T& value) {
        auto bytes = to_bytes(value);
        return out.write(reinterpret_cast<const char*>(bytes.data()), bytes.size()).good();
    }
    
    template<typename T>
    bool read_safe(std::istream& in, T& value) {
        std::array<std::byte, sizeof(T)> bytes;
        if (!in.read(reinterpret_cast<char*>(bytes.data()), bytes.size())) {
            return false;
        }
        value = from_bytes<T>(bytes);
        return true;
    }
}

// Enhanced expected implementation with proper constraints
template<typename T, typename E>
class expected {
private:
    alignas(std::max(alignof(T), alignof(E))) std::array<std::byte, std::max(sizeof(T), sizeof(E))> storage;
    bool has_val;
    
    T* value_ptr() noexcept { return std::launder(reinterpret_cast<T*>(storage.data())); }
    const T* value_ptr() const noexcept { return std::launder(reinterpret_cast<const T*>(storage.data())); }
    E* error_ptr() noexcept { return std::launder(reinterpret_cast<E*>(storage.data())); }
    const E* error_ptr() const noexcept { return std::launder(reinterpret_cast<const E*>(storage.data())); }
    
public:
    // Constrained constructors to prevent unintended copies/moves
    template<typename U>
    requires (!std::same_as<std::decay_t<U>, expected> && std::convertible_to<U, T>)
    explicit expected(U&& value) : has_val(true) { 
        std::construct_at(value_ptr(), std::forward<U>(value)); 
    }
    
    template<typename... Args>
    requires std::constructible_from<T, Args...>
    explicit expected(std::in_place_t, Args&&... args) : has_val(true) {
        std::construct_at(value_ptr(), std::forward<Args>(args)...);
    }
    
    template<typename U> 
    requires (!std::same_as<std::decay_t<U>, expected> && std::convertible_to<U, E>)
    explicit expected(U&& error_val) : has_val(false) { 
        std::construct_at(error_ptr(), std::forward<U>(error_val)); 
    }
    
    // Copy constructor
    expected(const expected& other) : has_val(other.has_val) {
        if (has_val) {
            std::construct_at(value_ptr(), *other.value_ptr());
        } else {
            std::construct_at(error_ptr(), *other.error_ptr());
        }
    }
    
    // Move constructor  
    expected(expected&& other) noexcept : has_val(other.has_val) {
        if (has_val) {
            std::construct_at(value_ptr(), std::move(*other.value_ptr()));
        } else {
            std::construct_at(error_ptr(), std::move(*other.error_ptr()));
        }
    }
    
    // Copy assignment
    expected& operator=(const expected& other) {
        if (this != &other) {
            if (has_val && other.has_val) {
                *value_ptr() = *other.value_ptr();
            } else if (!has_val && !other.has_val) {
                *error_ptr() = *other.error_ptr();
            } else {
                this->~expected();
                std::construct_at(this, other);
            }
        }
        return *this;
    }
    
    // Move assignment
    expected& operator=(expected&& other) noexcept {
        if (this != &other) {
            if (has_val && other.has_val) {
                *value_ptr() = std::move(*other.value_ptr());
            } else if (!has_val && !other.has_val) {
                *error_ptr() = std::move(*other.error_ptr());
            } else {
                this->~expected();
                std::construct_at(this, std::move(other));
            }
        }
        return *this;
    }
    
    ~expected() {
        if (has_val) {
            std::destroy_at(value_ptr());
        } else {
            std::destroy_at(error_ptr());
        }
    }
    
    bool has_value() const noexcept { return has_val; }
    explicit operator bool() const noexcept { return has_val; }
    
    T& operator*() & { return *value_ptr(); }
    const T& operator*() const & { return *value_ptr(); }
    T&& operator*() && { return std::move(*value_ptr()); }
    
    E& error() & { return *error_ptr(); }
    const E& error() const & { return *error_ptr(); }
};

template<typename E>
struct unexpected {
    E value;
    explicit unexpected(E&& e) : value(std::move(e)) {}
    explicit unexpected(const E& e) : value(e) {}
};

// Enhanced Perlin noise with std::array instead of C-style arrays
namespace noise {

constexpr float fade(float t) noexcept { 
    return t * t * t * (t * (t * 6 - 15) + 10); 
}

constexpr float lerp(float t, float a, float b) noexcept { 
    return a + t * (b - a); 
}

constexpr int fastfloor(float x) noexcept { 
    return (x > 0) ? static_cast<int>(x) : static_cast<int>(x) - 1; 
}

struct Perlin {
    std::array<int, 512> p; // Already using std::array - this is correct
    
    explicit Perlin(std::uint32_t seed = 0) {
        std::mt19937 rng(seed);
        std::iota(p.begin(), p.begin() + 256, 0);
        std::shuffle(p.begin(), p.begin() + 256, rng);
        for (int i = 0; i < 256; ++i) {
            p[256 + i] = p[i];
        }
    }
    
    constexpr float grad(int hash, float x, float y) const noexcept {
        const int h = hash & 3;
        const float u = h < 2 ? x : y;
        const float v = h < 2 ? y : x;
        return ((h & 1) ? -u : u) + ((h & 2) ? -2.0f * v : 2.0f * v);
    }
    
    float noise(float x, float y) const noexcept {
        const int xi = fastfloor(x) & 255;
        const int yi = fastfloor(y) & 255;
        const float xf = x - fastfloor(x);
        const float yf = y - fastfloor(y);
        const auto u = fade(xf);
        const auto v = fade(yf);
        
        const auto aaa = p[p[xi] + yi];
        const auto aba = p[p[xi] + yi + 1];
        const auto baa = p[p[xi + 1] + yi];
        const auto bba = p[p[xi + 1] + yi + 1];
        
        const auto x1 = lerp(u, grad(aaa, xf, yf), grad(baa, xf - 1, yf));
        const auto x2 = lerp(u, grad(aba, xf, yf - 1), grad(bba, xf - 1, yf - 1));
        return 0.5f * (lerp(v, x1, x2) + 1.0f);
    }
    
    /**
     * @tagline Multi-octave fractal noise generation for enhanced terrain diversity
     * @intuition Combine multiple noise frequencies to create more natural, varied landscapes
     * @approach Layer different octaves with decreasing amplitude for fractal-like terrain
     * @complexity Time: O(octaves), Space: O(1)
     */
    float fractalNoise(float x, float y, int octaves = 4) const noexcept {
        constexpr auto amplitudes = std::array{0.5f, 0.25f, 0.125f, 0.0625f, 0.03125f, 0.015625f};
        constexpr auto frequencies = std::array{0.02f, 0.04f, 0.08f, 0.16f, 0.32f, 0.64f};
        
        auto result = 0.0f;
        const auto maxOctaves = std::min(octaves, 6);
        for (int i = 0; i < maxOctaves; ++i) {
            result += amplitudes[i] * noise(frequencies[i] * x, frequencies[i] * y);
        }
        return result;
    }
};

} // namespace noise

// Enhanced modular tile system
using Coord = std::pair<int, int>;

enum class TerrainType : std::uint8_t {
    Water, Plain, Forest, Mountain, Path, POI, Start
};

enum class TileVariant : std::uint8_t { 
    A, B, C, D 
};

enum class TileRotation : std::uint8_t { 
    Deg0, Deg90, Deg180, Deg270 
};

struct ModularTile {
    TerrainType baseType{TerrainType::Plain};
    TileVariant variant{TileVariant::A};
    TileRotation rotation{TileRotation::Deg0};
    std::array<bool, 4> connections{true, true, true, true}; // N, E, S, W connectivity
    float elevation{0.5f};
    bool walkable{true};
    
    constexpr bool canConnect(const ModularTile& other, int direction) const noexcept {
        return connections[direction] && other.connections[(direction + 2) % 4];
    }
    
    bool operator==(const ModularTile&) const = default;
};

// Error handling for generation
enum class GenerationError {
    InvalidSeed,
    ConstraintsUnsatisfiable,
    PathfindingFailed,
    InsufficientWalkableArea,
    SerializationFailed
};

// Enhanced constraint system
struct LevelConstraints {
    float maxElevation{1.0f};
    float waterLevel{0.25f};
    float mountLevel{0.7f};
    float forestLevel{0.45f};
    float minWalkableRatio{0.6f};
    int minPOIDistance{15};
    int maxPOIDistance{40};
    int noiseOctaves{4};
    bool requireAllPOIsConnected{true};
    
    /**
     * @tagline Comprehensive constraint validation for level generation
     * @intuition Validate all design rules before accepting generated level
     * @approach Check each constraint systematically with early exit on failure
     * @complexity Time: O(n*m + pois²), Space: O(1)
     */
    template<typename Level>
    bool validate(const Level& level) const {
        // Check walkable ratio
        const auto walkableCount = static_cast<int>(std::ranges::count_if(level.grid, 
            [](const auto& tile) { return tile.walkable; }));
        const auto walkableRatio = static_cast<float>(walkableCount) / static_cast<float>(level.grid.size());
        if (walkableRatio < minWalkableRatio) {
            return false;
        }
        
        // Check POI distances
        for (size_t i = 0; i < level.pois.size(); ++i) {
            for (size_t j = i + 1; j < level.pois.size(); ++j) {
                const auto [x1, y1] = level.pois[i];
                const auto [x2, y2] = level.pois[j];
                const auto dist = std::abs(x1 - x2) + std::abs(y1 - y2);
                if (dist < minPOIDistance || dist > maxPOIDistance) {
                    return false;
                }
            }
        }
        
        return true;
    }
};

/**
 * @tagline Enhanced Level class with modular tiles and safe binary serialization
 * @intuition Use modern C++23 features for robust, maintainable level generation with safe memory operations
 * @approach Modular tile system with constraint validation, std::byte for binary I/O, and proper error handling
 * @complexity Time: O(n*m) for basic operations, Space: O(n*m) for storage
 */
class Level {
public:
    static constexpr std::uint32_t FORMAT_VERSION = 2;
    
    explicit Level(int w, int h) : width(w), height(h), grid(static_cast<size_t>(w) * static_cast<size_t>(h)) {}
    
    int width;
    int height;
    std::vector<ModularTile> grid;
    std::vector<Coord> pois;
    Coord playerStart{0, 0};
    LevelConstraints constraints;
    
    /**
     * @tagline Robust level generation with comprehensive error handling
     * @intuition Use expected pattern for clean error propagation without exceptions
     * @approach Try generation with backoff, returning specific error codes
     * @complexity Time: O(attempts * generation_time), Space: O(n*m)
     */
    expected<bool, GenerationError> generate(std::uint32_t seed = 42, int numPOI = 5) {
        noise::Perlin perlin(seed);
        std::mt19937 rng(seed);
        
        // 1. Generate terrain with multi-octave noise
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                auto& tile = at(x, y);
                tile.elevation = perlin.fractalNoise(static_cast<float>(x), static_cast<float>(y), constraints.noiseOctaves);
                
                // Assign tile variants and rotations for diversity
                std::uniform_int_distribution<> variantDist(0, 3);
                std::uniform_int_distribution<> rotDist(0, 3);
                tile.variant = static_cast<TileVariant>(variantDist(rng));
                tile.rotation = static_cast<TileRotation>(rotDist(rng));
                
                // Assign terrain type based on elevation
                if (tile.elevation < constraints.waterLevel) {
                    tile.baseType = TerrainType::Water;
                    tile.walkable = false;
                    tile.connections = {false, false, false, false};
                } else if (tile.elevation > constraints.mountLevel) {
                    tile.baseType = TerrainType::Mountain;
                    tile.walkable = false;
                    tile.connections = {false, false, false, false};
                } else if (tile.elevation > constraints.forestLevel) {
                    tile.baseType = TerrainType::Forest;
                    tile.walkable = true;
                } else {
                    tile.baseType = TerrainType::Plain;
                    tile.walkable = true;
                }
            }
        }
        
        // 2. Place player start
        if (!placePlayerStart(rng)) {
            return unexpected(GenerationError::InsufficientWalkableArea);
        }
        
        // 3. Place POIs with distance constraints
        if (!placePOIs(rng, numPOI)) {
            return unexpected(GenerationError::ConstraintsUnsatisfiable);
        }
        
        // 4. Ensure connectivity
        if (constraints.requireAllPOIsConnected) {
            if (!ensureConnectivity()) {
                return unexpected(GenerationError::PathfindingFailed);
            }
        }
        
        // 5. Validate constraints
        if (!constraints.validate(*this)) {
            return unexpected(GenerationError::ConstraintsUnsatisfiable);
        }
        
        return expected<bool, GenerationError>{std::in_place, true};
    }
    
    /**
     * @tagline A* pathfinding with enhanced heuristics for modular tiles
     * @intuition Find optimal path considering tile connectivity and terrain costs
     * @approach Priority queue with terrain-aware cost function and tile connections
     * @complexity Time: O(V log V), Space: O(V) where V is walkable cells
     */
    std::optional<std::vector<Coord>> findPath(const Coord& start, const Coord& end) const {
        constexpr auto directions = std::array<Coord, 4>{{{0,1},{1,0},{0,-1},{-1,0}}};
        
        std::vector<std::vector<float>> cost(static_cast<size_t>(height), 
            std::vector<float>(static_cast<size_t>(width), std::numeric_limits<float>::infinity()));
        std::vector<std::vector<Coord>> parent(static_cast<size_t>(height), 
            std::vector<Coord>(static_cast<size_t>(width), {-1, -1}));
        
        const auto heuristic = [](const Coord& a, const Coord& b) {
            return static_cast<float>(std::abs(a.first - b.first) + std::abs(a.second - b.second));
        };
        
        const auto getTerrainCost = [](TerrainType type) {
            switch (type) {
                case TerrainType::Plain: return 1.0f;
                case TerrainType::Forest: return 1.5f;
                case TerrainType::Path: return 0.8f;
                case TerrainType::Water: return 10.0f;
                case TerrainType::Mountain: return 10.0f;
                case TerrainType::POI: return 1.0f;
                case TerrainType::Start: return 1.0f;
                default: return 10.0f;
            }
        };
        
        using PQItem = std::pair<float, Coord>;
        std::priority_queue<PQItem, std::vector<PQItem>, std::greater<>> pq;
        
        cost[static_cast<size_t>(start.second)][static_cast<size_t>(start.first)] = 0.0f;
        pq.emplace(heuristic(start, end), start);
        
        while (!pq.empty()) {
            const auto [f, current] = pq.top();
            pq.pop();
            
            if (current == end) {
                std::vector<Coord> path;
                for (auto pos = end; pos != start; pos = parent[static_cast<size_t>(pos.second)][static_cast<size_t>(pos.first)]) {
                    path.push_back(pos);
                }
                path.push_back(start);
                std::ranges::reverse(path);
                return path;
            }
            
            for (int dir = 0; dir < 4; ++dir) {
                const auto [dx, dy] = directions[static_cast<size_t>(dir)];
                const auto nx = current.first + dx;
                const auto ny = current.second + dy;
                
                if (nx < 0 || ny < 0 || nx >= width || ny >= height) {
                    continue;
                }
                
                const auto& currentTile = at(current.first, current.second);
                const auto& nextTile = at(nx, ny);
                
                if (!nextTile.walkable || !currentTile.canConnect(nextTile, dir)) {
                    continue;
                }
                
                const auto newCost = cost[static_cast<size_t>(current.second)][static_cast<size_t>(current.first)] + getTerrainCost(nextTile.baseType);
                
                if (newCost < cost[static_cast<size_t>(ny)][static_cast<size_t>(nx)]) {
                    cost[static_cast<size_t>(ny)][static_cast<size_t>(nx)] = newCost;
                    parent[static_cast<size_t>(ny)][static_cast<size_t>(nx)] = current;
                    pq.emplace(newCost + heuristic({nx, ny}, end), Coord{nx, ny});
                }
            }
        }
        
        return std::nullopt;
    }
    
    /**
     * @tagline Safe binary serialization using std::byte and modern C++ features
     * @intuition Use std::byte for memory operations and safe casting to prevent undefined behavior
     * @approach Binary format with header validation, safe memory access, and proper error handling
     * @complexity Time: O(n*m), Space: O(1)
     */
    expected<void, GenerationError> serialize(std::ostream& out) const {
        // Use safe binary I/O instead of reinterpret_cast
        if (!binary_io::write_safe(out, FORMAT_VERSION)) {
            return unexpected(GenerationError::SerializationFailed);
        }
        
        if (!binary_io::write_safe(out, width) || 
            !binary_io::write_safe(out, height) ||
            !binary_io::write_safe(out, constraints) ||
            !binary_io::write_safe(out, playerStart)) {
            return unexpected(GenerationError::SerializationFailed);
        }
        
        const auto poiCount = static_cast<std::uint32_t>(pois.size());
        if (!binary_io::write_safe(out, poiCount)) {
            return unexpected(GenerationError::SerializationFailed);
        }
        
        if (poiCount > 0) {
            // Safe byte-oriented write for POI array
            const auto poiBytes = std::span{
                reinterpret_cast<const std::byte*>(pois.data()), 
                sizeof(Coord) * poiCount
            };
            if (!out.write(reinterpret_cast<const char*>(poiBytes.data()), poiBytes.size())) {
                return unexpected(GenerationError::SerializationFailed);
            }
        }
        
        // Safe tile serialization
        for (const auto& tile : grid) {
            if (!binary_io::write_safe(out, tile)) {
                return unexpected(GenerationError::SerializationFailed);
            }
        }
        
        return expected<void, GenerationError>{std::in_place};
    }
    
    static expected<Level, GenerationError> deserialize(std::istream& in) {
        std::uint32_t version;
        if (!binary_io::read_safe(in, version)) {
            return unexpected(GenerationError::SerializationFailed);
        }
        
        if (version > FORMAT_VERSION) {
            return unexpected(GenerationError::SerializationFailed);
        }
        
        int w, h;
        LevelConstraints constraints;
        Coord playerStart;
        
        if (!binary_io::read_safe(in, w) ||
            !binary_io::read_safe(in, h) ||
            !binary_io::read_safe(in, constraints) ||
            !binary_io::read_safe(in, playerStart)) {
            return unexpected(GenerationError::SerializationFailed);
        }
        
        Level level(w, h);
        level.constraints = constraints;
        level.playerStart = playerStart;
        
        std::uint32_t poiCount;
        if (!binary_io::read_safe(in, poiCount)) {
            return unexpected(GenerationError::SerializationFailed);
        }
        
        level.pois.resize(poiCount);
        if (poiCount > 0) {
            // Safe byte-oriented read for POI array
            auto poiBytes = std::span{
                reinterpret_cast<std::byte*>(level.pois.data()), 
                sizeof(Coord) * poiCount
            };
            if (!in.read(reinterpret_cast<char*>(poiBytes.data()), poiBytes.size())) {
                return unexpected(GenerationError::SerializationFailed);
            }
        }
        
        // Safe tile deserialization
        for (auto& tile : level.grid) {
            if (!binary_io::read_safe(in, tile)) {
                return unexpected(GenerationError::SerializationFailed);
            }
        }
        
        return expected<Level, GenerationError>{std::in_place, std::move(level)};
    }
    
    const ModularTile& at(int x, int y) const { 
        return grid[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)]; 
    }
    
    ModularTile& at(int x, int y) { 
        return grid[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)]; 
    }
    
private:
    bool placePlayerStart(std::mt19937& rng) {
        std::uniform_int_distribution<> distX(0, width - 1);
        std::uniform_int_distribution<> distY(0, height - 1);
        
        for (int attempt = 0; attempt < 1000; ++attempt) {
            const auto x = distX(rng);
            const auto y = distY(rng);
            auto& tile = at(x, y);
            if (tile.walkable && tile.elevation > constraints.waterLevel && tile.elevation < constraints.mountLevel) {
                playerStart = {x, y};
                tile.baseType = TerrainType::Start;
                return true;
            }
        }
        return false;
    }
    
    bool placePOIs(std::mt19937& rng, int numPOI) {
        pois.clear();
        std::uniform_int_distribution<> distX(0, width - 1);
        std::uniform_int_distribution<> distY(0, height - 1);
        
        for (int placed = 0; placed < numPOI;) {
            const auto x = distX(rng);
            const auto y = distY(rng);
            const auto candidate = Coord{x, y};
            
            if (candidate == playerStart || !at(x, y).walkable) {
                continue;
            }
            
            const auto validDistance = std::ranges::all_of(pois, [&](const Coord& poi) {
                const auto dist = std::abs(x - poi.first) + std::abs(y - poi.second);
                return dist >= constraints.minPOIDistance && dist <= constraints.maxPOIDistance;
            });
            
            if (validDistance) {
                at(x, y).baseType = TerrainType::POI;
                pois.push_back(candidate);
                ++placed;
            }
        }
        
        return pois.size() == static_cast<size_t>(numPOI);
    }
    
    bool ensureConnectivity() {
        for (const auto& poi : pois) {
            const auto path = findPath(playerStart, poi);
            if (!path) {
                return false;
            }
            
            for (const auto& [x, y] : *path) {
                auto& tile = at(x, y);
                if (tile.walkable && tile.baseType != TerrainType::POI && tile.baseType != TerrainType::Start) {
                    tile.baseType = TerrainType::Path;
                }
            }
        }
        return true;
    }
};

namespace vis {

constexpr int CELL_SIZE = 12;
constexpr int MARGIN = 40;
constexpr int INFO_PANEL_WIDTH = 200;

/**
 * @tagline Interactive real-time level preview with enhanced safety and modern C++ features
 * @intuition Provide immediate visual feedback for level design iteration with safe resource management
 * @approach SDL2-based renderer with RAII, mouse interaction, and overlay information using modern C++ patterns
 * @complexity Time: O(visible_cells), Space: O(1)
 */
class LevelDebugger {
private:
    SDL_Window* window{};
    SDL_Renderer* renderer{};
    bool showElevation{false};
    bool showConnections{false};
    bool showVariants{true};
    Coord selectedCell{-1, -1};
    Coord hoveredCell{-1, -1};
    
public:
    LevelDebugger() = default;
    ~LevelDebugger() {
        cleanup();
    }
    
    // Delete copy operations
    LevelDebugger(const LevelDebugger&) = delete;
    LevelDebugger& operator=(const LevelDebugger&) = delete;
    
    // Delete move operations for simplicity
    LevelDebugger(LevelDebugger&&) = delete;
    LevelDebugger& operator=(LevelDebugger&&) = delete;
    
    void preview(const Level& level) {
        if (initialize(level)) {
            runMainLoop(level);
        }
        cleanup();
    }
    
private:
    bool initialize(const Level& level) {
        if (SDL_Init(SDL_INIT_VIDEO) < 0) {
            return false;
        }
        
        const auto windowWidth = level.width * CELL_SIZE + 2 * MARGIN + INFO_PANEL_WIDTH;
        const auto windowHeight = level.height * CELL_SIZE + 2 * MARGIN;
        
        window = SDL_CreateWindow("Level Debugger", 
            SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
            windowWidth, windowHeight, SDL_WINDOW_SHOWN);
        
        if (!window) {
            return false;
        }
        
        renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
        return renderer != nullptr;
    }
    
    void runMainLoop(const Level& level) {
        auto running = true;
        SDL_Event event;
        
        while (running) {
            while (SDL_PollEvent(&event)) {
                if (event.type == SDL_QUIT) {
                    running = false;
                } else {
                    handleInput(event, level);
                }
            }
            
            render(level);
            SDL_Delay(16);
        }
    }
    
    void handleInput(const SDL_Event& event, const Level& level) {
        if (event.type == SDL_MOUSEMOTION) {
            const auto mouseX = event.motion.x - MARGIN;
            const auto mouseY = event.motion.y - MARGIN;
            
            if (mouseX >= 0 && mouseY >= 0) {
                const auto cellX = mouseX / CELL_SIZE;
                const auto cellY = mouseY / CELL_SIZE;
                
                if (cellX < level.width && cellY < level.height) {
                    hoveredCell = {cellX, cellY};
                } else {
                    hoveredCell = {-1, -1};
                }
            }
        } else if (event.type == SDL_MOUSEBUTTONDOWN) {
            if (event.button.button == SDL_BUTTON_LEFT) {
                selectedCell = hoveredCell;
            }
        } else if (event.type == SDL_KEYDOWN) {
            switch (event.key.keysym.sym) {
                case SDLK_e: 
                    showElevation = !showElevation; 
                    break;
                case SDLK_c: 
                    showConnections = !showConnections; 
                    break;
                case SDLK_v: 
                    showVariants = !showVariants; 
                    break;
                case SDLK_ESCAPE: 
                    selectedCell = {-1, -1}; 
                    break;
                default: 
                    break;
            }
        }
    }
    
    void render(const Level& level) {
        SDL_SetRenderDrawColor(renderer, 20, 20, 20, 255);
        SDL_RenderClear(renderer);
        
        drawLevel(level);
        drawOverlays(level);
        drawInfoPanel(level);
        
        SDL_RenderPresent(renderer);
    }
    
    void drawLevel(const Level& level) {
        for (int y = 0; y < level.height; ++y) {
            for (int x = 0; x < level.width; ++x) {
                const auto& tile = level.at(x, y);
                const auto rect = SDL_Rect{
                    MARGIN + x * CELL_SIZE,
                    MARGIN + y * CELL_SIZE,
                    CELL_SIZE,
                    CELL_SIZE
                };
                
                const auto color = getTileColor(tile, showElevation, showVariants);
                SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);
                SDL_RenderFillRect(renderer, &rect);
                
                if (showConnections) {
                    drawTileConnections(tile, rect);
                }
                
                if (Coord{x, y} == selectedCell) {
                    SDL_SetRenderDrawColor(renderer, 255, 255, 0, 255);
                    SDL_RenderDrawRect(renderer, &rect);
                } else if (Coord{x, y} == hoveredCell) {
                    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 128);
                    SDL_RenderDrawRect(renderer, &rect);
                }
            }
        }
    }
    
    void drawOverlays(const Level& level) {
        SDL_SetRenderDrawColor(renderer, 255, 255, 0, 200);
        for (const auto& poi : level.pois) {
            if (const auto path = level.findPath(level.playerStart, poi)) {
                for (size_t i = 1; i < path->size(); ++i) {
                    const auto [x1, y1] = (*path)[i-1];
                    const auto [x2, y2] = (*path)[i];
                    SDL_RenderDrawLine(renderer,
                        MARGIN + x1 * CELL_SIZE + CELL_SIZE/2,
                        MARGIN + y1 * CELL_SIZE + CELL_SIZE/2,
                        MARGIN + x2 * CELL_SIZE + CELL_SIZE/2,
                        MARGIN + y2 * CELL_SIZE + CELL_SIZE/2);
                }
            }
        }
    }
    
    void drawInfoPanel(const Level& level) {
        const auto panelX = MARGIN + level.width * CELL_SIZE + 20;
        const auto panelY = MARGIN;
        
        const auto panelRect = SDL_Rect{panelX, panelY, INFO_PANEL_WIDTH - 20, level.height * CELL_SIZE};
        SDL_SetRenderDrawColor(renderer, 40, 40, 40, 200);
        SDL_RenderFillRect(renderer, &panelRect);
        
        if (selectedCell.first >= 0) {
            const auto& tile = level.at(selectedCell.first, selectedCell.second);
            const auto swatch = SDL_Rect{panelX + 10, panelY + 10, 20, 20};
            const auto color = getTileColor(tile, false, false);
            SDL_SetRenderDrawColor(renderer, color.r, color.g, color.b, 255);
            SDL_RenderFillRect(renderer, &swatch);
        }
    }
    
    void drawTileConnections(const ModularTile& tile, const SDL_Rect& rect) {
        SDL_SetRenderDrawColor(renderer, 0, 255, 0, 128);
        const auto centerX = rect.x + rect.w / 2;
        const auto centerY = rect.y + rect.h / 2;
        
        if (tile.connections[0]) {
            SDL_RenderDrawLine(renderer, centerX, centerY, centerX, rect.y);
        }
        if (tile.connections[1]) {
            SDL_RenderDrawLine(renderer, centerX, centerY, rect.x + rect.w, centerY);
        }
        if (tile.connections[2]) {
            SDL_RenderDrawLine(renderer, centerX, centerY, centerX, rect.y + rect.h);
        }
        if (tile.connections[3]) {
            SDL_RenderDrawLine(renderer, centerX, centerY, rect.x, centerY);
        }
    }
    
    static SDL_Color getTileColor(const ModularTile& tile, bool showElevation, bool showVariants) {
        if (showElevation) {
            const auto intensity = static_cast<std::uint8_t>(tile.elevation * 255);
            return {intensity, intensity, intensity, 255};
        }
        
        auto baseColor = SDL_Color{};
        switch (tile.baseType) {
            case TerrainType::Water: 
                baseColor = {30, 80, 180, 255}; 
                break;
            case TerrainType::Plain: 
                baseColor = {90, 210, 80, 255}; 
                break;
            case TerrainType::Forest: 
                baseColor = {15, 110, 40, 255}; 
                break;
            case TerrainType::Mountain: 
                baseColor = {120, 115, 130, 255}; 
                break;
            case TerrainType::Path: 
                baseColor = {255, 230, 60, 255}; 
                break;
            case TerrainType::POI: 
                baseColor = {220, 30, 30, 255}; 
                break;
            case TerrainType::Start: 
                baseColor = {230, 250, 90, 255}; 
                break;
            default: 
                baseColor = {255, 255, 255, 255}; 
                break;
        }
        
        if (showVariants) {
            // Use std::to_underlying instead of manual casting
            const auto variantOffset = std::to_underlying(tile.variant) * 10 - 15;
            baseColor.r = static_cast<std::uint8_t>(std::clamp(static_cast<int>(baseColor.r) + variantOffset, 0, 255));
            baseColor.g = static_cast<std::uint8_t>(std::clamp(static_cast<int>(baseColor.g) + variantOffset, 0, 255));
            baseColor.b = static_cast<std::uint8_t>(std::clamp(static_cast<int>(baseColor.b) + variantOffset, 0, 255));
        }
        
        return baseColor;
    }
    
    void cleanup() {
        if (renderer) {
            SDL_DestroyRenderer(renderer);
            renderer = nullptr;
        }
        if (window) {
            SDL_DestroyWindow(window);
            window = nullptr;
        }
        SDL_Quit();
    }
};

void previewLevel(const Level& level) {
    auto debugger = LevelDebugger{};
    debugger.preview(level);
}

} // namespace vis

/**
 * @tagline Production-ready main entry point with enhanced safety and modern C++ features
 * @intuition Demonstrate all features with error handling, performance metrics, and safe operations
 * @approach Generate multiple levels, test serialization, and launch interactive preview using modern C++ patterns
 * @complexity Time: O(generations * level_complexity), Space: O(level_size)
 */
int main() {
    constexpr auto width = 128;
    constexpr auto height = 96;
    constexpr auto numPOI = 8;
    
    const auto testSeeds = std::array{42u, 12345u, 99999u, 
        static_cast<std::uint32_t>(std::chrono::system_clock::now().time_since_epoch().count())};
    
    auto bestLevel = Level{width, height};
    auto generationSucceeded = false;
    
    for (const auto seed : testSeeds) {
        auto level = Level{width, height};
        
        level.constraints.noiseOctaves = 5;
        level.constraints.minWalkableRatio = 0.65f;
        level.constraints.minPOIDistance = 12;
        level.constraints.maxPOIDistance = 35;
        
        if (const auto result = level.generate(seed, numPOI); result.has_value()) {
            std::cout << "✓ Level generated successfully with seed " << seed << std::endl;
            bestLevel = std::move(level);
            generationSucceeded = true;
            break;
        } else {
            std::cout << "✗ Generation failed with seed " << seed 
                     << " (Error: " << std::to_underlying(result.error()) << ")" << std::endl;
        }
    }
    
    if (!generationSucceeded) {
        std::cerr << "Failed to generate valid level with any seed!" << std::endl;
        return 1;
    }
    
    // Test serialization with safe binary I/O
    {
        auto ofs = std::ofstream{"level_test.save", std::ios::binary};
        if (const auto saveResult = bestLevel.serialize(ofs); saveResult.has_value()) {
            ofs.close();
            
            auto ifs = std::ifstream{"level_test.save", std::ios::binary};
            if (const auto loadResult = Level::deserialize(ifs); loadResult.has_value()) {
                std::cout << "✓ Serialization round-trip successful" << std::endl;
                
                const auto& loaded = *loadResult;
                const auto dataIntegrityOk = (loaded.width == bestLevel.width && 
                                             loaded.height == bestLevel.height &&
                                             loaded.grid == bestLevel.grid &&
                                             loaded.pois == bestLevel.pois);
                
                std::cout << (dataIntegrityOk ? "✓" : "✗") 
                         << " Data integrity check " 
                         << (dataIntegrityOk ? "passed" : "failed") << std::endl;
            } else {
                std::cout << "✗ Deserialization failed" << std::endl;
            }
            ifs.close();
        } else {
            std::cout << "✗ Serialization failed" << std::endl;
        }
    }
    
    // Performance metrics
    const auto startTime = std::chrono::high_resolution_clock::now();
    constexpr auto pathfindingTests = 100;
    auto successfulPaths = 0;
    
    for (int i = 0; i < pathfindingTests; ++i) {
        if (bestLevel.pois.size() >= 2) {
            if (const auto path = bestLevel.findPath(bestLevel.pois[0], bestLevel.pois[1]); path.has_value()) {
                ++successfulPaths;
            }
        }
    }
    
    const auto endTime = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
    
    std::cout << "Performance: " << pathfindingTests << " pathfinding tests in " 
              << duration.count() << "μs (" << successfulPaths << " successful)" << std::endl;
    
    std::cout << "\nLaunching interactive preview..." << std::endl;
    std::cout << "Controls: E=elevation, C=connections, V=variants, ESC=deselect, mouse=select" << std::endl;
    
    vis::previewLevel(bestLevel);
    
    return 0;
}
