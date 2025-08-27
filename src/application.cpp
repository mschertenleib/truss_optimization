#include "application.hpp"
#include "analysis.hpp"
#include "unique_resource.hpp"
#include "vec.hpp"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define GLFW_INCLUDE_ES3
#define GL_GLES_PROTOTYPES 0
#else
#define GLFW_INCLUDE_GLCOREARB
#endif
#define GLFW_INCLUDE_GLEXT
#include <GLFW/glfw3.h>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <optional>
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
    f(PFNGLGETUNIFORMLOCATIONPROC, glGetUniformLocation);                      \
    f(PFNGLUNIFORM2FPROC, glUniform2f);                                        \
    f(PFNGLGENBUFFERSPROC, glGenBuffers);                                      \
    f(PFNGLDELETEBUFFERSPROC, glDeleteBuffers);                                \
    f(PFNGLBINDBUFFERPROC, glBindBuffer);                                      \
    f(PFNGLBUFFERDATAPROC, glBufferData);                                      \
    f(PFNGLBUFFERSUBDATAPROC, glBufferSubData);                                \
    f(PFNGLGENVERTEXARRAYSPROC, glGenVertexArrays);                            \
    f(PFNGLDELETEVERTEXARRAYSPROC, glDeleteVertexArrays);                      \
    f(PFNGLBINDVERTEXARRAYPROC, glBindVertexArray);                            \
    f(PFNGLUNIFORMBLOCKBINDINGPROC, glUniformBlockBinding);                    \
    f(PFNGLVERTEXATTRIBPOINTERPROC, glVertexAttribPointer);                    \
    f(PFNGLENABLEVERTEXATTRIBARRAYPROC, glEnableVertexAttribArray);            \
    f(PFNGLDRAWELEMENTSPROC, glDrawElements);                                  \
    f(PFNGLUSEPROGRAMPROC, glUseProgram);                                      \
    f(PFNGLVIEWPORTPROC, glViewport);                                          \
    f(PFNGLCLEARCOLORPROC, glClearColor);                                      \
    f(PFNGLCLEARPROC, glClear);                                                \
    f(PFNGLBLENDFUNCPROC, glBlendFunc);

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

struct ImGui_deleter
{
    void operator()(Stateless)
    {
        ImGui::DestroyContext();
    }
};

struct ImGui_glfw_deleter
{
    void operator()(Stateless)
    {
        ImGui_ImplGlfw_Shutdown();
    }
};

struct ImGui_opengl_deleter
{
    void operator()(Stateless)
    {
        ImGui_ImplOpenGL3_Shutdown();
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

struct Window_state
{
    float scale_x;
    float scale_y;
    int framebuffer_width;
    int framebuffer_height;
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
};

struct Circle
{
    vec2 center;
    float radius;
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

struct Application_state
{
    Analysis_state analysis;
    Unique_resource<Stateless, GLFW_deleter> glfw_context;
    Unique_resource<GLFWwindow *, Window_deleter> window;
    Unique_resource<Stateless, ImGui_deleter> imgui_context;
    Unique_resource<Stateless, ImGui_glfw_deleter> imgui_glfw_context;
    Unique_resource<Stateless, ImGui_opengl_deleter> imgui_opengl_context;
    Window_state window_state;
    std::vector<Line> lines;
    std::vector<Circle> circles;
    Raster_geometry raster_geometry;
    Unique_resource<GLuint, GL_array_deleter> vao {};
    Unique_resource<GLuint, GL_array_deleter> vbo {};
    Unique_resource<GLuint, GL_array_deleter> ibo {};
    Unique_resource<GLuint, GL_deleter> line_program {};
    Unique_resource<GLuint, GL_deleter> circle_program {};
    GLint loc_view_position_draw_line {};
    GLint loc_view_size_draw_line {};
    GLint loc_view_position_draw_circle {};
    GLint loc_view_size_draw_circle {};
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

void glfw_error_callback(int error, const char *description)
{
    std::cerr << "GLFW error " << error << ": " << description << '\n';
}

void glfw_window_content_scale_callback(GLFWwindow *window,
                                        float xscale,
                                        float yscale)
{
    auto *const window_state =
        static_cast<Window_state *>(glfwGetWindowUserPointer(window));
    assert(window_state != nullptr);

    window_state->scale_x = xscale;
    window_state->scale_y = yscale;
}

void glfw_framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
    auto *const window_state =
        static_cast<Window_state *>(glfwGetWindowUserPointer(window));
    assert(window_state != nullptr);

    window_state->framebuffer_width = width;
    window_state->framebuffer_height = height;
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
    shader_code.push_back("\n");
    if (std::string_view(glsl_version).ends_with("es"))
    {
        shader_code.push_back("precision highp float;\n");
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

uniform vec2 view_position;
uniform vec2 view_size;

out vec4 local;
out vec3 color;

void main()
{
    vec2 position = (vertex_position - view_position) / (0.5 * view_size);
    gl_Position = vec4(position, 0.0, 1.0);
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

void create_raster_geometry(const std::vector<Line> &lines,
                            const std::vector<Circle> &circles,
                            float thickness,
                            Raster_geometry &geometry)
{
    geometry.vertices.clear();
    geometry.indices.clear();

    constexpr vec3 color {1.0f, 1.0f, 1.0f};

    geometry.line_indices_offset = geometry.indices.size();
    for (const auto &line : lines)
    {
        const auto line_vec = line.b - line.a;
        const auto line_length = norm(line_vec);
        const auto line_dir = line_vec * (1.0f / line_length);
        const auto delta_left =
            vec2 {-line_dir.y, line_dir.x} * (thickness * 0.5f);
        const auto delta_up = line_dir * (thickness * 0.5f);
        const auto start_left = line.a + delta_left - delta_up;
        const auto start_right = line.a - delta_left - delta_up;
        const auto end_left = line.b + delta_left + delta_up;
        const auto end_right = line.b - delta_left + delta_up;
        const auto aspect_ratio = line_length / thickness;

        const auto first_index =
            static_cast<std::uint32_t>(geometry.vertices.size());
        geometry.vertices.push_back(
            {start_left, {-0.5f, 0.5f, -aspect_ratio - 0.5f, 0.0f}, color});
        geometry.vertices.push_back(
            {start_right, {0.5f, 0.5f, -aspect_ratio - 0.5f, 0.0f}, color});
        geometry.vertices.push_back(
            {end_right, {0.5f, -aspect_ratio - 0.5f, 0.5f, 0.0f}, color});
        geometry.vertices.push_back(
            {end_left, {-0.5f, -aspect_ratio - 0.5f, 0.5f, 0.0f}, color});
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
        const auto half_side = circle.radius + 0.5f * thickness;
        const auto bottom_left = circle.center + vec2 {-half_side, -half_side};
        const auto bottom_right = circle.center + vec2 {half_side, -half_side};
        const auto top_right = circle.center + vec2 {half_side, half_side};
        const auto top_left = circle.center + vec2 {-half_side, half_side};
        const auto rel_thickness = thickness / half_side;

        const auto first_index =
            static_cast<std::uint32_t>(geometry.vertices.size());
        geometry.vertices.push_back(
            {bottom_left, {-1.0, -1.0f, rel_thickness, 0.0f}, color});
        geometry.vertices.push_back(
            {bottom_right, {1.0f, -1.0f, rel_thickness, 0.0f}, color});
        geometry.vertices.push_back(
            {top_right, {1.0f, 1.0f, rel_thickness, 0.0f}, color});
        geometry.vertices.push_back(
            {top_left, {-1.0f, 1.0f, rel_thickness, 0.0f}, color});
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

[[nodiscard]] Application_state init_application()
{
    Application_state app {};

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

    // FIXME: title
    auto *const window_ptr =
        glfwCreateWindow(1280, 720, "OLC CodeJam 2025", nullptr, nullptr);
    if (window_ptr == nullptr)
    {
        throw std::runtime_error("Failed to create GLFW window");
    }
    app.window.reset(window_ptr);

    glfwMakeContextCurrent(app.window.get());

    glfwSetWindowUserPointer(app.window.get(), &app.window_state);
    glfwSetWindowContentScaleCallback(app.window.get(),
                                      &glfw_window_content_scale_callback);
    glfwSetFramebufferSizeCallback(app.window.get(),
                                   &glfw_framebuffer_size_callback);
    glfwGetWindowContentScale(
        app.window.get(), &app.window_state.scale_x, &app.window_state.scale_y);
    glfwGetFramebufferSize(app.window.get(),
                           &app.window_state.framebuffer_width,
                           &app.window_state.framebuffer_height);

    glfwSwapInterval(1);

    load_gl_functions();

#ifndef __EMSCRIPTEN__
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
    glDebugMessageCallback(&gl_debug_callback, nullptr);
#endif

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    app.imgui_context.reset(Stateless {});

    ImGui::GetIO().ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    ImGui::StyleColorsDark();

    if (!ImGui_ImplGlfw_InitForOpenGL(app.window.get(), true))
    {
        throw std::runtime_error("ImGui: failed to initialize GLFW backend");
    }
    app.imgui_glfw_context.reset(Stateless {});

#ifdef __EMSCRIPTEN__
    ImGui_ImplGlfw_InstallEmscriptenCallbacks(app.window.get(), "#canvas");
#endif

    if (!ImGui_ImplOpenGL3_Init(glsl_version))
    {
        throw std::runtime_error("ImGui: failed to initialize OpenGL backend");
    }
    app.imgui_opengl_context.reset(Stateless {});

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Generate truss
    constexpr std::uint32_t num_nodes {21};
    float x {-1.2f};
    float y {0.0f};
    constexpr float dx {2.4f / (num_nodes - 1)};
    constexpr float dy {0.4f};
    for (std::uint32_t i {0}; i < num_nodes; ++i)
    {
        app.analysis.nodes.push_back({x, y});
        x += dx;
        y += i % 2 == 0 ? dy : -dy;
    }
    for (std::uint32_t i {0}; i < num_nodes - 1; ++i)
    {
        app.analysis.elements.push_back({i, i + 1});
        if (i < num_nodes - 2)
        {
            app.analysis.elements.push_back({i, i + 2});
        }
    }
    app.analysis.fixed_dofs = {0, 1, 2, 3};
    for (std::uint32_t i {4}; i < num_nodes * 2; ++i)
    {
        app.analysis.free_dofs.push_back(i);
    }
    app.analysis.loads.setZero(app.analysis.free_dofs.size());
    app.analysis.loads(num_nodes * 2 - 1 - 4) = -1.0f;

    app.lines.clear();
    for (const auto [i, j] : app.analysis.elements)
    {
        app.lines.push_back({app.analysis.nodes[i], app.analysis.nodes[j]});
    }
    create_raster_geometry(app.lines, app.circles, 0.03f, app.raster_geometry);

    // TODO: we should probably organize the raster geometry better. We need to
    // clarify the distinction between updating the vertex buffer (when moving
    // objects) and re-creating the vertex and index buffers with a new size
    // (when adding or removing objects).
    std::tie(app.vao, app.vbo, app.ibo) =
        create_vertex_index_buffers(app.raster_geometry);

    app.line_program = create_program(
        glsl_version, vertex_shader_code, line_fragment_shader_code);
    app.circle_program = create_program(
        glsl_version, vertex_shader_code, circle_fragment_shader_code);
    app.loc_view_position_draw_line =
        glGetUniformLocation(app.line_program.get(), "view_position");
    app.loc_view_size_draw_line =
        glGetUniformLocation(app.line_program.get(), "view_size");
    app.loc_view_position_draw_circle =
        glGetUniformLocation(app.circle_program.get(), "view_position");
    app.loc_view_size_draw_circle =
        glGetUniformLocation(app.circle_program.get(), "view_size");

    return app;
}

void draw_geometry(const Application_state &app)
{
    // FIXME
    const float view_x {0.0f};
    const float view_y {0.0f};
    const float view_height {2.0f};
    const float view_width {
        view_height * static_cast<float>(app.window_state.framebuffer_width) /
        static_cast<float>(app.window_state.framebuffer_height)};

    glBindVertexArray(app.vao.get());

    glUseProgram(app.line_program.get());
    glUniform2f(app.loc_view_position_draw_line, view_x, view_y);
    glUniform2f(app.loc_view_size_draw_line, view_width, view_height);
    glDrawElements(
        GL_TRIANGLES,
        static_cast<GLsizei>(app.raster_geometry.line_indices_size),
        GL_UNSIGNED_INT,
        reinterpret_cast<void *>(app.raster_geometry.line_indices_offset *
                                 sizeof(std::uint32_t)));

    glUseProgram(app.circle_program.get());
    glUniform2f(app.loc_view_position_draw_circle, view_x, view_y);
    glUniform2f(app.loc_view_size_draw_circle, view_width, view_height);
    glDrawElements(
        GL_TRIANGLES,
        static_cast<GLsizei>(app.raster_geometry.circle_indices_size),
        GL_UNSIGNED_INT,
        reinterpret_cast<void *>(app.raster_geometry.circle_indices_offset *
                                 sizeof(std::uint32_t)));

    glBindVertexArray(0);
}

void main_loop_update(Application_state &app)
{
    glfwPollEvents();

    static float force_factor {};
    const auto base_loads = app.analysis.loads;
    app.analysis.loads *= force_factor;
    assemble(app.analysis);
    solve(app.analysis);
    app.analysis.loads = base_loads;

    app.lines.clear();
    for (const auto [i, j] : app.analysis.elements)
    {
        const auto node_i = app.analysis.nodes[i];
        const vec2 disp_i {app.analysis.displacements(2 * i),
                           app.analysis.displacements(2 * i + 1)};
        const auto node_j = app.analysis.nodes[j];
        const vec2 disp_j {app.analysis.displacements(2 * j),
                           app.analysis.displacements(2 * j + 1)};
        app.lines.push_back({node_i + disp_i, node_j + disp_j});
    }
    create_raster_geometry(app.lines, app.circles, 0.03f, app.raster_geometry);
    update_vertex_buffer(app.vao.get(), app.vbo.get(), app.raster_geometry);

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    if (ImGui::Begin("Test"))
    {
        ImGui::Text("%.2f ms/frame, %.2f fps",
                    static_cast<double>(1000.0f / ImGui::GetIO().Framerate),
                    static_cast<double>(ImGui::GetIO().Framerate));
        ImGui::SliderFloat("Force", &force_factor, 0.0f, 500000.0f);
    }
    ImGui::End();

    ImGui::Render();

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    draw_geometry(app);

    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(app.window.get());
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

    ImGui::GetIO().IniFilename = nullptr;

    emscripten_set_main_loop_arg(
        [](void *arg)
        { main_loop_update(*static_cast<Application_state *>(arg)); },
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
