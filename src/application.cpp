#include "application.hpp"
#include "optimization.hpp"
#include "unique_resource.hpp"
#include "vec.hpp"

#define STB_TRUETYPE_IMPLEMENTATION
#define STBTT_STATIC
#include "stb_truetype.h"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define GLFW_INCLUDE_ES3
#define GL_GLES_PROTOTYPES 0
#include <GLFW/emscripten_glfw3.h>
#else
#define GLFW_INCLUDE_GLCOREARB
#endif
#define GLFW_INCLUDE_GLEXT
#include <GLFW/glfw3.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

#define ENUMERATE_GL_FUNCTIONS_COMMON(f)                                       \
    f(PFNGLENABLEPROC, glEnable);                                              \
    f(PFNGLCREATESHADERPROC, glCreateShader);                                  \
    f(PFNGLDELETESHADERPROC, glDeleteShader);                                  \
    f(PFNGLSHADERSOURCEPROC, glShaderSource);                                  \
    f(PFNGLCOMPILESHADERPROC, glCompileShader);                                \
    f(PFNGLGETSHADERIVPROC, glGetShaderiv);                                    \
    f(PFNGLGETSHADERINFOLOGPROC, glGetShaderInfoLog);                          \
    f(PFNGLATTACHSHADERPROC, glAttachShader);                                  \
    f(PFNGLCREATEPROGRAMPROC, glCreateProgram);                                \
    f(PFNGLDELETEPROGRAMPROC, glDeleteProgram);                                \
    f(PFNGLLINKPROGRAMPROC, glLinkProgram);                                    \
    f(PFNGLGETPROGRAMIVPROC, glGetProgramiv);                                  \
    f(PFNGLGETPROGRAMINFOLOGPROC, glGetProgramInfoLog);                        \
    f(PFNGLGENBUFFERSPROC, glGenBuffers);                                      \
    f(PFNGLDELETEBUFFERSPROC, glDeleteBuffers);                                \
    f(PFNGLBINDBUFFERPROC, glBindBuffer);                                      \
    f(PFNGLBUFFERDATAPROC, glBufferData);                                      \
    f(PFNGLBUFFERSUBDATAPROC, glBufferSubData);                                \
    f(PFNGLGENVERTEXARRAYSPROC, glGenVertexArrays);                            \
    f(PFNGLDELETEVERTEXARRAYSPROC, glDeleteVertexArrays);                      \
    f(PFNGLBINDVERTEXARRAYPROC, glBindVertexArray);                            \
    f(PFNGLVERTEXATTRIBPOINTERPROC, glVertexAttribPointer);                    \
    f(PFNGLENABLEVERTEXATTRIBARRAYPROC, glEnableVertexAttribArray);            \
    f(PFNGLDRAWELEMENTSPROC, glDrawElements);                                  \
    f(PFNGLUSEPROGRAMPROC, glUseProgram);                                      \
    f(PFNGLVIEWPORTPROC, glViewport);                                          \
    f(PFNGLCLEARCOLORPROC, glClearColor);                                      \
    f(PFNGLCLEARPROC, glClear);                                                \
    f(PFNGLBLENDFUNCPROC, glBlendFunc);                                        \
    f(PFNGLBLENDEQUATIONPROC, glBlendEquation);                                \
    f(PFNGLGENTEXTURESPROC, glGenTextures);                                    \
    f(PFNGLDELETETEXTURESPROC, glDeleteTextures);                              \
    f(PFNGLBINDTEXTUREPROC, glBindTexture);                                    \
    f(PFNGLPIXELSTOREIPROC, glPixelStorei);                                    \
    f(PFNGLTEXIMAGE2DPROC, glTexImage2D);                                      \
    f(PFNGLTEXPARAMETERIPROC, glTexParameteri);

#define ENUMERATE_GL_FUNCTIONS_430(f)                                          \
    f(PFNGLDEBUGMESSAGECALLBACKPROC, glDebugMessageCallback);

#ifndef __EMSCRIPTEN__
#define ENUMERATE_GL_FUNCTIONS(f)                                              \
    ENUMERATE_GL_FUNCTIONS_COMMON(f) ENUMERATE_GL_FUNCTIONS_430(f)
#else
#define ENUMERATE_GL_FUNCTIONS(f) ENUMERATE_GL_FUNCTIONS_COMMON(f)
#endif

// clang-format off
#define DECLARE_GL_FUNCTION(type, name) type name {nullptr}
// clang-format on

ENUMERATE_GL_FUNCTIONS(DECLARE_GL_FUNCTION)

struct Stateless
{
};

struct GLFW_deleter
{
    void operator()(Stateless)
    {
        glfwTerminate();
    }
};

struct Window_deleter
{
    void operator()(GLFWwindow *window)
    {
        glfwDestroyWindow(window);
    }
};

struct GL_deleter
{
    void (*destroy)(GLuint);
    void operator()(GLuint handle)
    {
        destroy(handle);
    }
};

struct GL_array_deleter
{
    void (*destroy)(GLsizei, const GLuint *);
    void operator()(GLuint handle)
    {
        destroy(1, &handle);
    }
};

struct Vertex
{
    vec2 position;
    vec4 local;
    vec3 color;
};

struct Line
{
    vec2 a;
    vec2 b;
    float thickness;
    vec3 color;
};

struct Circle
{
    vec2 center;
    float radius;
    float thickness;
    vec3 color;
};

struct Raster_geometry
{
    std::size_t line_indices_offset;
    std::size_t line_indices_size;
    std::size_t circle_indices_offset;
    std::size_t circle_indices_size;
    std::vector<Vertex> vertices;
    std::vector<std::uint32_t> indices;
};

struct Viewport
{
    int x;
    int y;
    int width;
    int height;
};

struct Rectangle
{
    vec2 pos;
    vec2 size;
};

enum struct State
{
    idle,
    placing_support,
    placing_load,
    running
};

struct Application
{
    Analysis_state analysis;
    unsigned int step;
    State state;
    bool should_idle;
    Unique_resource<Stateless, GLFW_deleter> glfw_context;
    Unique_resource<GLFWwindow *, Window_deleter> window;
    float scale_x;
    float scale_y;
    int framebuffer_width;
    int framebuffer_height;
    Viewport viewport;
    std::vector<Line> lines;
    std::vector<Circle> circles;
    Raster_geometry raster_geometry;
    Unique_resource<GLuint, GL_array_deleter> vao {};
    Unique_resource<GLuint, GL_array_deleter> vbo {};
    Unique_resource<GLuint, GL_array_deleter> ibo {};
    Unique_resource<GLuint, GL_deleter> line_program {};
    Unique_resource<GLuint, GL_deleter> circle_program {};
    static constexpr vec2 world_center {0.0f, 0.0f};
    static constexpr vec2 world_size {2.0f, 2.0f};
    static constexpr Rectangle play_button {0.85f, 0.15f, 0.1f, 0.1f};
    static constexpr Rectangle step_button {0.85f, 0.0f, 0.1f, 0.1f};
    static constexpr Rectangle restart_button {0.85f, -0.15f, 0.1f, 0.1f};
};

template <std::invocable C, std::invocable<GLuint> D>
[[nodiscard]] auto create_object(C &&create, D &&destroy)
{
    return Unique_resource(create(), GL_deleter {destroy});
}

template <std::invocable<GLenum> C, std::invocable<GLuint> D>
[[nodiscard]] auto create_object(C &&create, GLenum arg, D &&destroy)
{
    return Unique_resource(create(arg), GL_deleter {destroy});
}

template <std::invocable<GLsizei, GLuint *> C,
          std::invocable<GLsizei, const GLuint *> D>
[[nodiscard]] auto create_object(C &&create, D &&destroy)
{
    GLuint object {};
    create(1, &object);
    return Unique_resource(object, GL_array_deleter {destroy});
}

[[nodiscard]] constexpr float screen_to_world(float x,
                                              int screen_min,
                                              int screen_size,
                                              float world_center,
                                              float world_size) noexcept
{
    const auto u =
        (x - static_cast<float>(screen_min)) / static_cast<float>(screen_size);
    return world_center + (u - 0.5f) * world_size;
}

[[nodiscard]] constexpr Viewport centered_viewport(
    float aspect_ratio, int available_width, int available_height) noexcept
{
    const auto available_aspect_ratio = static_cast<float>(available_width) /
                                        static_cast<float>(available_height);
    if (aspect_ratio > available_aspect_ratio)
    {
        const auto height = static_cast<int>(
            static_cast<float>(available_width) / aspect_ratio);
        return {0, (available_height - height) / 2, available_width, height};
    }
    else
    {
        const auto width = static_cast<int>(
            static_cast<float>(available_height) * aspect_ratio);
        return {(available_width - width) / 2, 0, width, available_height};
    }
}

[[nodiscard]] constexpr bool
is_in_rectangle(const vec2 &point, const Rectangle &rectangle) noexcept
{
    return point.x >= rectangle.pos.x &&
           point.x <= rectangle.pos.x + rectangle.size.x &&
           point.y >= rectangle.pos.y &&
           point.y <= rectangle.pos.y + rectangle.size.y;
}

void glfw_error_callback(int error, const char *description)
{
    std::cerr << "GLFW error " << error << ": " << description << '\n';
}

void glfw_mouse_button_callback(GLFWwindow *window,
                                int button,
                                int action,
                                [[maybe_unused]] int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        auto *const app =
            static_cast<Application *>(glfwGetWindowUserPointer(window));
        assert(app != nullptr);

        double mouse_x {};
        double mouse_y {};
        glfwGetCursorPos(app->window.get(), &mouse_x, &mouse_y);
        const auto mouse_world_x =
            screen_to_world(static_cast<float>(mouse_x) * app->scale_x,
                            app->viewport.x,
                            app->viewport.width,
                            app->world_center.x,
                            app->world_size.x);
        const auto mouse_world_y =
            screen_to_world(static_cast<float>(mouse_y) * app->scale_y,
                            app->viewport.y + app->viewport.height,
                            -app->viewport.height,
                            app->world_center.y,
                            app->world_size.y);
        const vec2 mouse_pos {mouse_world_x, mouse_world_y};

        if (is_in_rectangle(mouse_pos, app->play_button))
        {
            if (app->state == State::idle)
            {
                app->state = State::running;
            }
            else if (app->state == State::running)
            {
                app->state = State::idle;
            }
        }
        else if (app->state == State::idle &&
                 is_in_rectangle(mouse_pos, app->step_button))
        {
            app->state = State::running;
            app->should_idle = true;
        }
        else if (is_in_rectangle(mouse_pos, app->restart_button))
        {
            app->step = 0;
            const std::vector<vec2> fixed_nodes {{-0.8f, 0.4f}, {-0.8f, -0.4f}};
            const vec2 load_node {0.8f, -0.2f};
            const vec2 load_vector {0.0f, -1.0f};
            optimization_init(
                fixed_nodes, load_node, load_vector, app->analysis);
        }
    }
}

void glfw_window_content_scale_callback(GLFWwindow *window,
                                        float xscale,
                                        float yscale)
{
    auto *const app =
        static_cast<Application *>(glfwGetWindowUserPointer(window));
    assert(app != nullptr);

    app->scale_x = xscale;
    app->scale_y = yscale;
}

void glfw_framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
    auto *const app =
        static_cast<Application *>(glfwGetWindowUserPointer(window));
    assert(app != nullptr);

    app->framebuffer_width = width;
    app->framebuffer_height = height;
    app->viewport = centered_viewport(1.0f, width, height);
}

void load_gl_functions()
{
#define LOAD_GL_FUNCTION(type, name)                                           \
    name = reinterpret_cast<type>(glfwGetProcAddress(#name));                  \
    assert(name != nullptr)

    ENUMERATE_GL_FUNCTIONS(LOAD_GL_FUNCTION)

#undef LOAD_GL_FUNCTION
}

#ifndef __EMSCRIPTEN__
void APIENTRY gl_debug_callback([[maybe_unused]] GLenum source,
                                GLenum type,
                                [[maybe_unused]] GLuint id,
                                GLenum severity,
                                [[maybe_unused]] GLsizei length,
                                const GLchar *message,
                                [[maybe_unused]] const void *user_param)
{
    if (type == GL_DEBUG_TYPE_OTHER ||
        severity == GL_DEBUG_SEVERITY_NOTIFICATION)
    {
        return;
    }
    std::cerr << message << '\n';
}
#endif

[[nodiscard]] auto
create_shader(GLenum type, std::size_t size, const char *const code[])
{
    auto shader = create_object(glCreateShader, type, glDeleteShader);

    glShaderSource(shader.get(), static_cast<GLsizei>(size), code, nullptr);
    glCompileShader(shader.get());

    int success {};
    glGetShaderiv(shader.get(), GL_COMPILE_STATUS, &success);
    if (!success)
    {
        int buf_length {};
        glGetShaderiv(shader.get(), GL_INFO_LOG_LENGTH, &buf_length);
        std::string message(static_cast<std::size_t>(buf_length), '\0');
        glGetShaderInfoLog(shader.get(), buf_length, nullptr, message.data());
        std::ostringstream oss;
        oss << "Shader compilation failed:\n" << message << '\n';
        throw std::runtime_error(oss.str());
    }

    return shader;
}

[[nodiscard]] auto create_program(const char *glsl_version,
                                  const char *vertex_shader_code,
                                  const char *fragment_shader_code)
{
    std::vector<const char *> shader_code;
    shader_code.reserve(3);
    shader_code.push_back(glsl_version);
    if (std::string_view(glsl_version).ends_with("es"))
    {
        shader_code.push_back("\nprecision highp float;\n");
    }
    else
    {
        shader_code.push_back("\n");
    }
    shader_code.push_back(vertex_shader_code);

    const auto vertex_shader =
        create_shader(GL_VERTEX_SHADER, shader_code.size(), shader_code.data());

    shader_code.back() = fragment_shader_code;
    const auto fragment_shader = create_shader(
        GL_FRAGMENT_SHADER, shader_code.size(), shader_code.data());

    auto program = create_object(glCreateProgram, glDeleteProgram);
    glAttachShader(program.get(), vertex_shader.get());
    glAttachShader(program.get(), fragment_shader.get());
    glLinkProgram(program.get());

    int success {};
    glGetProgramiv(program.get(), GL_LINK_STATUS, &success);
    if (!success)
    {
        int buf_length {};
        glGetProgramiv(program.get(), GL_INFO_LOG_LENGTH, &buf_length);
        std::string message(static_cast<std::size_t>(buf_length), '\0');
        glGetProgramInfoLog(program.get(), buf_length, nullptr, message.data());
        std::ostringstream oss;
        oss << "Program linking failed:\n" << message << '\n';
        throw std::runtime_error(oss.str());
    }

    return program;
}

constexpr auto vertex_shader_code = R"(
layout (location = 0) in vec2 vertex_position;
layout (location = 1) in vec4 vertex_local;
layout (location = 2) in vec3 vertex_color;

out vec4 local;
out vec3 color;

void main()
{
    gl_Position = vec4(vertex_position, 0.0, 1.0);
    local = vertex_local;
    color = vertex_color;
})";

constexpr auto line_fragment_shader_code = R"(
in vec4 local;
in vec3 color;

out vec4 frag_color;

void main()
{
    float dist = abs(local.x);
    if (local.y > 0.0 || local.z > 0.0)
    {
        dist = min(length(local.xy), length(local.xz));
    }
    // NOTE: this assumes that the local x, y and z have the same scale
    float pixel_size = fwidth(local.x);
    float alpha = clamp((0.5 - dist) / pixel_size, 0.0, 1.0);
    frag_color = vec4(color, alpha);
})";

constexpr auto circle_fragment_shader_code = R"(
in vec4 local;
in vec3 color;

out vec4 frag_color;

void main()
{
    float dist = length(local.xy);
    float thickness = local.z;
    float pixel_size = fwidth(dist);
    float alpha = clamp((1.0 - dist) / pixel_size, 0.0, 1.0)
        - clamp((1.0 - thickness + pixel_size - dist) / pixel_size, 0.0, 1.0);
    frag_color = vec4(color, alpha);
})";

[[nodiscard]] auto create_vertex_index_buffers(const Raster_geometry &geometry)
{
    auto vao = create_object(glGenVertexArrays, glDeleteVertexArrays);
    glBindVertexArray(vao.get());

    auto vbo = create_object(glGenBuffers, glDeleteBuffers);
    glBindBuffer(GL_ARRAY_BUFFER, vbo.get());
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizei>(geometry.vertices.size() * sizeof(Vertex)),
        geometry.vertices.data(),
        GL_DYNAMIC_DRAW);

    auto ibo = create_object(glGenBuffers, glDeleteBuffers);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo.get());
    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER,
        static_cast<GLsizei>(geometry.indices.size() * sizeof(std::uint32_t)),
        geometry.indices.data(),
        GL_STATIC_DRAW);

    glVertexAttribPointer(0,
                          sizeof(Vertex::position) / sizeof(float),
                          GL_FLOAT,
                          GL_FALSE,
                          sizeof(Vertex),
                          reinterpret_cast<void *>(offsetof(Vertex, position)));
    glVertexAttribPointer(1,
                          sizeof(Vertex::local) / sizeof(float),
                          GL_FLOAT,
                          GL_FALSE,
                          sizeof(Vertex),
                          reinterpret_cast<void *>(offsetof(Vertex, local)));
    glVertexAttribPointer(2,
                          sizeof(Vertex::color) / sizeof(float),
                          GL_FLOAT,
                          GL_FALSE,
                          sizeof(Vertex),
                          reinterpret_cast<void *>(offsetof(Vertex, color)));
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);

    glBindVertexArray(0);

    return std::tuple {std::move(vao), std::move(vbo), std::move(ibo)};
}

void update_vertex_buffer(GLuint vao,
                          GLuint vbo,
                          const Raster_geometry &geometry)
{
    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferSubData(
        GL_ARRAY_BUFFER,
        0,
        static_cast<GLsizei>(geometry.vertices.size() * sizeof(Vertex)),
        geometry.vertices.data());

    glBindVertexArray(0);
}

[[nodiscard]] std::vector<std::uint8_t> read_binary_file(const char *file_name)
{
    const std::filesystem::path path(file_name);
    if (!std::filesystem::exists(path))
    {
        std::ostringstream oss;
        oss << "File " << path << " does not exist";
        throw std::runtime_error(oss.str());
    }

    const auto file_size = std::filesystem::file_size(path);
    std::vector<std::uint8_t> buffer(file_size);

    std::ifstream file(path, std::ios::binary);
    if (!file)
    {
        std::ostringstream oss;
        oss << "Failed to open file " << path;
        throw std::runtime_error(oss.str());
    }

    file.read(reinterpret_cast<char *>(buffer.data()),
              static_cast<std::streamsize>(file_size));
    if (file.eof())
    {
        std::ostringstream oss;
        oss << "End-of-file reached while reading file " << path;
        throw std::runtime_error(oss.str());
    }

    return buffer;
}

void create_font_texture()
{
    const auto ttf = read_binary_file("font/Roboto_Condensed-Regular.ttf");
    constexpr int atlas_width {1024};
    constexpr int atlas_height {1024};
    std::vector<std::uint8_t> pixels(atlas_width * atlas_height);

    stbtt_pack_context pc {};
    stbtt_PackBegin(
        &pc, pixels.data(), atlas_width, atlas_height, 0, 1, nullptr);
    stbtt_PackSetOversampling(&pc, 2, 2);
    constexpr int num_chars {95};
    stbtt_packedchar char_data[num_chars] {};
    stbtt_PackFontRange(&pc, ttf.data(), 0, 32.0f, 32, num_chars, char_data);
    stbtt_PackEnd(&pc);

    auto font_texture = create_object(glGenTextures, glDeleteTextures);
    glBindTexture(GL_TEXTURE_2D, font_texture.get());
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_R8,
                 atlas_width,
                 atlas_height,
                 0,
                 GL_RED,
                 GL_UNSIGNED_BYTE,
                 pixels.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    /*
    // --- When drawing a string ---
float x = 20.0f, y = 50.0f; // baseline position in screen space
stbtt_aligned_quad q;
for (const char* p = "Hello OpenGL!"; *p; ++p) {
    if (*p < 32 || *p > 126) continue;
    stbtt_GetPackedQuad(cdata, ATLAS_W, ATLAS_H, *p - 32, &x, &y, &q, 1);
    // append 4 verts (q.x0/q.y0..q.x1/q.y1, q.s0/q.t0..q.s1/q.t1) to your
vertex buffer
}
// In your fragment shader: float a = texture(uAtlas, uv).r; outColor =
vec4(textColor, a);
// Make sure blending is enabled: SRC_ALPHA, ONE_MINUS_SRC_ALPHA.
    */
}

void create_raster_geometry(const std::vector<Line> &lines,
                            const std::vector<Circle> &circles,
                            Raster_geometry &geometry)
{
    geometry.vertices.clear();
    geometry.indices.clear();

    geometry.line_indices_offset = geometry.indices.size();
    for (const auto &line : lines)
    {
        const auto line_vec = line.b - line.a;
        const auto line_length = norm(line_vec);
        const auto line_dir = line_vec * (1.0f / line_length);
        const auto delta_left =
            vec2 {-line_dir.y, line_dir.x} * (line.thickness * 0.5f);
        const auto delta_up = line_dir * (line.thickness * 0.5f);
        const auto start_left = line.a + delta_left - delta_up;
        const auto start_right = line.a - delta_left - delta_up;
        const auto end_left = line.b + delta_left + delta_up;
        const auto end_right = line.b - delta_left + delta_up;
        const auto aspect_ratio = line_length / line.thickness;

        const auto first_index =
            static_cast<std::uint32_t>(geometry.vertices.size());
        geometry.vertices.push_back({start_left,
                                     {-0.5f, 0.5f, -aspect_ratio - 0.5f, 0.0f},
                                     line.color});
        geometry.vertices.push_back({start_right,
                                     {0.5f, 0.5f, -aspect_ratio - 0.5f, 0.0f},
                                     line.color});
        geometry.vertices.push_back(
            {end_right, {0.5f, -aspect_ratio - 0.5f, 0.5f, 0.0f}, line.color});
        geometry.vertices.push_back(
            {end_left, {-0.5f, -aspect_ratio - 0.5f, 0.5f, 0.0f}, line.color});
        geometry.indices.push_back(first_index + 0);
        geometry.indices.push_back(first_index + 1);
        geometry.indices.push_back(first_index + 2);
        geometry.indices.push_back(first_index + 0);
        geometry.indices.push_back(first_index + 2);
        geometry.indices.push_back(first_index + 3);
    }
    geometry.line_indices_size =
        geometry.indices.size() - geometry.line_indices_offset;

    geometry.circle_indices_offset = geometry.indices.size();
    for (const auto &circle : circles)
    {
        const auto half_side = circle.radius + 0.5f * circle.thickness;
        const auto bottom_left = circle.center + vec2 {-half_side, -half_side};
        const auto bottom_right = circle.center + vec2 {half_side, -half_side};
        const auto top_right = circle.center + vec2 {half_side, half_side};
        const auto top_left = circle.center + vec2 {-half_side, half_side};
        const auto rel_thickness = circle.thickness / half_side;

        const auto first_index =
            static_cast<std::uint32_t>(geometry.vertices.size());
        geometry.vertices.push_back(
            {bottom_left, {-1.0, -1.0f, rel_thickness, 0.0f}, circle.color});
        geometry.vertices.push_back(
            {bottom_right, {1.0f, -1.0f, rel_thickness, 0.0f}, circle.color});
        geometry.vertices.push_back(
            {top_right, {1.0f, 1.0f, rel_thickness, 0.0f}, circle.color});
        geometry.vertices.push_back(
            {top_left, {-1.0f, 1.0f, rel_thickness, 0.0f}, circle.color});
        geometry.indices.push_back(first_index + 0);
        geometry.indices.push_back(first_index + 1);
        geometry.indices.push_back(first_index + 2);
        geometry.indices.push_back(first_index + 0);
        geometry.indices.push_back(first_index + 2);
        geometry.indices.push_back(first_index + 3);
    }
    geometry.circle_indices_size =
        geometry.indices.size() - geometry.circle_indices_offset;
}

[[nodiscard]] Application init_application()
{
    Application app {};

    glfwSetErrorCallback(&glfw_error_callback);

    if (!glfwInit())
    {
        throw std::runtime_error("Failed to initialize GLFW");
    }
    app.glfw_context.reset(Stateless {});

#ifdef __EMSCRIPTEN__
    // WebGL 2.0
    constexpr auto glsl_version = "#version 300 es";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_ES_API);
#else
    constexpr auto glsl_version = "#version 430 core";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_CONTEXT_DEBUG, GLFW_TRUE);
#endif
    glfwWindowHint(GLFW_SAMPLES, 0);

    auto *const window_ptr =
        glfwCreateWindow(1280, 720, "Truss Optimization", nullptr, nullptr);
    if (window_ptr == nullptr)
    {
        throw std::runtime_error("Failed to create GLFW window");
    }
    app.window.reset(window_ptr);

    glfwMakeContextCurrent(app.window.get());

    glfwSetWindowUserPointer(app.window.get(), &app);

    glfwSetMouseButtonCallback(app.window.get(), &glfw_mouse_button_callback);
    glfwSetWindowContentScaleCallback(app.window.get(),
                                      &glfw_window_content_scale_callback);
    glfwSetFramebufferSizeCallback(app.window.get(),
                                   &glfw_framebuffer_size_callback);

    glfwGetWindowContentScale(app.window.get(), &app.scale_x, &app.scale_y);
    glfwGetFramebufferSize(
        app.window.get(), &app.framebuffer_width, &app.framebuffer_height);
    app.viewport =
        centered_viewport(1.0f, app.framebuffer_width, app.framebuffer_height);

    glfwSwapInterval(1);

    load_gl_functions();

#ifndef __EMSCRIPTEN__
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
    glDebugMessageCallback(&gl_debug_callback, nullptr);
#endif

#ifdef __EMSCRIPTEN__
    emscripten_glfw_make_canvas_resizable(app.window.get(), "window", nullptr);
#endif

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    app.line_program = create_program(
        glsl_version, vertex_shader_code, line_fragment_shader_code);
    app.circle_program = create_program(
        glsl_version, vertex_shader_code, circle_fragment_shader_code);

    // create_font_texture();

    app.state = State::idle;

    const std::vector<vec2> fixed_nodes {{-0.8f, 0.4f}, {-0.8f, -0.4f}};
    const vec2 load_node {0.8f, -0.2f};
    const vec2 load_vector {0.0f, -1.0f};
    optimization_init(fixed_nodes, load_node, load_vector, app.analysis);

    return app;
}

void draw_geometry(const Application &app)
{
    glBindVertexArray(app.vao.get());

    glUseProgram(app.line_program.get());
    glDrawElements(
        GL_TRIANGLES,
        static_cast<GLsizei>(app.raster_geometry.line_indices_size),
        GL_UNSIGNED_INT,
        reinterpret_cast<void *>(app.raster_geometry.line_indices_offset *
                                 sizeof(std::uint32_t)));

    glUseProgram(app.circle_program.get());
    glDrawElements(
        GL_TRIANGLES,
        static_cast<GLsizei>(app.raster_geometry.circle_indices_size),
        GL_UNSIGNED_INT,
        reinterpret_cast<void *>(app.raster_geometry.circle_indices_offset *
                                 sizeof(std::uint32_t)));

    glBindVertexArray(0);
}

void main_loop_update(Application &app)
{
    glfwPollEvents();

    if (app.state == State::running)
    {
        optimization_step(app.analysis);
        ++app.step;

        if (app.should_idle)
        {
            app.state = State::idle;
            app.should_idle = false;
        }
    }

    const auto min_force = app.analysis.axial_forces.minCoeff();
    const auto max_force = app.analysis.axial_forces.maxCoeff();
    app.lines.clear();
    for (std::size_t e {0}; e < app.analysis.elements.size(); ++e)
    {
        const auto [i, j] = app.analysis.elements[e];
        const auto activation =
            app.analysis.activations[static_cast<Eigen::Index>(e)];
        const auto force =
            app.analysis.axial_forces[static_cast<Eigen::Index>(e)];

        const auto rel_force =
            force >= 0.0f ? force / max_force : force / min_force;
        const auto max_color = force >= 0.0f ? vec3 {0.25f, 0.25f, 1.0f}
                                             : vec3 {1.0f, 0.25f, 0.25f};
        const auto color = rel_force * max_color + 1.0f - rel_force;
        app.lines.push_back({.a = app.analysis.nodes[i],
                             .b = app.analysis.nodes[j],
                             .thickness = activation * 0.03f,
                             .color = color});
    }

    constexpr vec3 color {0.75f, 0.75f, 0.75f};
    constexpr float thickness {0.02f};

    for (const auto &rectangle :
         {app.play_button, app.step_button, app.restart_button})
    {
        const auto p0 = rectangle.pos;
        const auto p1 = rectangle.pos + vec2 {rectangle.size.x, 0.0f};
        const auto p2 = rectangle.pos + rectangle.size;
        const auto p3 = rectangle.pos + vec2 {0.0f, rectangle.size.y};
        app.lines.emplace_back(p0, p1, thickness, color);
        app.lines.emplace_back(p1, p2, thickness, color);
        app.lines.emplace_back(p2, p3, thickness, color);
        app.lines.emplace_back(p3, p0, thickness, color);
    }

    if (app.state == State::running)
    {
        const auto p0 =
            app.play_button.pos + app.play_button.size * vec2 {0.3f, 0.3f};
        const auto p1 =
            app.play_button.pos + app.play_button.size * vec2 {0.3f, 0.7f};
        const auto p2 =
            app.play_button.pos + app.play_button.size * vec2 {0.7f, 0.3f};
        const auto p3 =
            app.play_button.pos + app.play_button.size * vec2 {0.7f, 0.7f};
        app.lines.emplace_back(p0, p1, thickness, color);
        app.lines.emplace_back(p2, p3, thickness, color);
    }
    else
    {
        const auto p0 =
            app.play_button.pos + app.play_button.size * vec2 {0.7f, 0.5f};
        const auto p1 =
            app.play_button.pos + app.play_button.size * vec2 {0.3f, 0.8f};
        const auto p2 =
            app.play_button.pos + app.play_button.size * vec2 {0.3f, 0.2f};
        app.lines.emplace_back(p0, p1, thickness, color);
        app.lines.emplace_back(p1, p2, thickness, color);
        app.lines.emplace_back(p2, p0, thickness, color);
    }

    const auto p0 =
        app.step_button.pos + app.step_button.size * vec2 {0.65f, 0.5f};
    const auto p1 =
        app.step_button.pos + app.step_button.size * vec2 {0.3f, 0.8f};
    const auto p2 =
        app.step_button.pos + app.step_button.size * vec2 {0.3f, 0.2f};
    const auto p3 =
        app.step_button.pos + app.step_button.size * vec2 {0.7f, 0.2f};
    const auto p4 =
        app.step_button.pos + app.step_button.size * vec2 {0.7f, 0.8f};
    app.lines.emplace_back(p0, p1, thickness, color);
    app.lines.emplace_back(p1, p2, thickness, color);
    app.lines.emplace_back(p2, p0, thickness, color);
    app.lines.emplace_back(p3, p4, thickness, color);

    app.circles.clear();
    for (std::size_t i {0}; i < app.analysis.nodes.size(); ++i)
    {
        app.circles.push_back({.center = app.analysis.nodes[i],
                               .radius = 0.02f,
                               .thickness = 0.02f,
                               .color = {1.0f, 1.0f, 1.0f}});
    }

    create_raster_geometry(app.lines, app.circles, app.raster_geometry);
    std::tie(app.vao, app.vbo, app.ibo) =
        create_vertex_index_buffers(app.raster_geometry);

    glViewport(app.viewport.x,
               app.viewport.y,
               app.viewport.width,
               app.viewport.height);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    draw_geometry(app);

#ifndef __EMSCRIPTEN__
    glfwSwapBuffers(app.window.get());
#endif

    const double now {glfwGetTime()};
    static double last_time {now};
    const auto elapsed = now - last_time;
    if (elapsed > 0.0)
    {
        std::cout << std::format("{:7.2f} ms, {:7.2f} fps\r",
                                 elapsed * 1000.0,
                                 1.0 / elapsed)
                  << std::flush;
        last_time = now;
    }
}

} // namespace

void run_application()
{
#ifdef __EMSCRIPTEN__

    // NOTE: The emscripten_set_main_loop(..., simulate_infinite_loop=true) call
    // unwinds the stack (but global or static storage duration objects are not
    // destroyed), before transferring control to the browser. For this reason,
    // any object passed to the main loop callback must have static storage
    // duration.

    static auto app = init_application();

    emscripten_set_main_loop_arg(
        [](void *arg) { main_loop_update(*static_cast<Application *>(arg)); },
        &app,
        0,
        true);

#else

    auto app = init_application();

    while (!glfwWindowShouldClose(app.window.get()))
    {
        main_loop_update(app);
    }

#endif
}
