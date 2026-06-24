#include "application.hpp"
#include "optimization.hpp"
#include "unique_handle.hpp"
#include "vec.hpp"

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

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
#include <format>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <optional>
#include <random>
#include <ranges>
#include <source_location>
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

#ifdef __EMSCRIPTEN__
#define ENUMERATE_GL_FUNCTIONS(f) ENUMERATE_GL_FUNCTIONS_COMMON(f)
#else
#define ENUMERATE_GL_FUNCTIONS(f)                                              \
    ENUMERATE_GL_FUNCTIONS_COMMON(f) ENUMERATE_GL_FUNCTIONS_430(f)
#endif

// clang-format off
#define DECLARE_GL_FUNCTION(type, name) type name {nullptr}
// clang-format on

ENUMERATE_GL_FUNCTIONS(DECLARE_GL_FUNCTION)

struct GLFW_deleter
{
    void operator()()
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
    void operator()(ImGuiContext *context)
    {
        ImGui::DestroyContext(context);
    }
};

struct ImGui_glfw_deleter
{
    void operator()()
    {
        ImGui_ImplGlfw_Shutdown();
    }
};

struct ImGui_opengl_deleter
{
    void operator()()
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

struct Viewport
{
    int x;
    int y;
    int width;
    int height;
};

enum struct App_state
{
    idle,
    placing_support,
    placing_load,
    running
};

struct Render_data
{
    std::vector<Line> lines;
    std::vector<Circle> circles;
    std::size_t line_indices_offset;
    std::size_t line_indices_size;
    std::size_t circle_indices_offset;
    std::size_t circle_indices_size;
    std::vector<Vertex> vertices;
    std::vector<std::uint32_t> indices;
    Unique_handle<GLuint, GL_array_deleter> vao;
    Unique_handle<GLuint, GL_array_deleter> vbo;
    Unique_handle<GLuint, GL_array_deleter> ibo;
    Unique_handle<GLuint, GL_deleter> line_program;
    Unique_handle<GLuint, GL_deleter> circle_program;
};

struct Application
{
    Optimization_state optimization;
    Problem problem;
    unsigned int step;
    App_state state;
    bool should_idle;
    Unique_handle<bool, GLFW_deleter> glfw_context;
    Unique_handle<GLFWwindow *, Window_deleter> window;
    Unique_handle<ImGuiContext *, ImGui_deleter> imgui_context;
    Unique_handle<bool, ImGui_glfw_deleter> imgui_glfw_context;
    Unique_handle<bool, ImGui_opengl_deleter> imgui_opengl_context;
    float scale_x;
    float scale_y;
    int framebuffer_width;
    int framebuffer_height;
    Viewport viewport;
    Render_data render_data;

    static constexpr vec2 world_center {0.0f, 0.0f};
    static constexpr vec2 world_size {2.0f, 2.0f};
};

inline void
check(bool condition,
      const char *expression,
      const std::source_location &loc = std::source_location::current())
{
    if (!condition)
    {
        throw std::runtime_error(std::format("{}:{}: check failed:\n  {}",
                                             loc.file_name(),
                                             loc.line(),
                                             expression));
    }
}

#define CHECK(expression) check(expression, #expression)

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

void glfw_error_callback(int error, const char *description)
{
    std::cerr << "GLFW error " << error << ": " << description << '\n';
}

void glfw_mouse_button_callback(GLFWwindow *window,
                                int button,
                                int action,
                                [[maybe_unused]] int mods)
{
    if (ImGui::GetIO().WantCaptureMouse)
    {
        return;
    }

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        auto *const app =
            static_cast<Application *>(glfwGetWindowUserPointer(window));
        assert(app != nullptr);

        double mouse_x {};
        double mouse_y {};
        glfwGetCursorPos(app->window.get(), &mouse_x, &mouse_y);
        const auto mouse_world_x = screen_to_world(static_cast<float>(mouse_x),
                                                   app->viewport.x,
                                                   app->viewport.width,
                                                   app->world_center.x,
                                                   app->world_size.x);
        const auto mouse_world_y =
            screen_to_world(static_cast<float>(mouse_y),
                            app->viewport.y + app->viewport.height,
                            -app->viewport.height,
                            app->world_center.y,
                            app->world_size.y);
        const vec2 mouse_pos {mouse_world_x, mouse_world_y};
        std::cout << std::format("Click [{}, {}]\n", mouse_pos.x, mouse_pos.y);
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
    app->viewport = centered_viewport(
        Application::world_size.x / Application::world_size.y, width, height);
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
    Unique_handle shader(glCreateShader(type), GL_deleter {glDeleteShader});

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
        throw std::runtime_error(
            std::format("Shader compilation failed:\n{}", message));
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
    const auto constants =
        std::format("const vec2 world_center = vec2({}, {});\n"
                    "const vec2 world_size = vec2({}, {});\n",
                    Application::world_center.x,
                    Application::world_center.y,
                    Application::world_size.x,
                    Application::world_size.y);
    shader_code.push_back(constants.c_str());
    shader_code.push_back(vertex_shader_code);

    const auto vertex_shader =
        create_shader(GL_VERTEX_SHADER, shader_code.size(), shader_code.data());

    shader_code.pop_back();
    shader_code.pop_back();
    shader_code.push_back(fragment_shader_code);
    const auto fragment_shader = create_shader(
        GL_FRAGMENT_SHADER, shader_code.size(), shader_code.data());

    Unique_handle program(glCreateProgram(), GL_deleter {glDeleteProgram});
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
        throw std::runtime_error(
            std::format("Program linking failed:\n{}", message));
    }

    return program;
}

inline constexpr auto vertex_shader_code = R"(
layout (location = 0) in vec2 vertex_position;
layout (location = 1) in vec4 vertex_local;
layout (location = 2) in vec3 vertex_color;

out vec4 local;
out vec3 color;

void main()
{
    vec2 position = (vertex_position - world_center) / (0.5 * world_size);
    gl_Position = vec4(position, 0.0, 1.0);
    local = vertex_local;
    color = vertex_color;
})";

inline constexpr auto line_fragment_shader_code = R"(
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

inline constexpr auto circle_fragment_shader_code = R"(
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

void create_vertex_and_index_buffers(Render_data &render_data)
{
    GLuint vao_gl {};
    glGenVertexArrays(1, &vao_gl);
    render_data.vao.reset(vao_gl, GL_array_deleter {glDeleteVertexArrays});
    glBindVertexArray(render_data.vao.get());

    GLuint vbo_gl {};
    glGenBuffers(1, &vbo_gl);
    render_data.vbo.reset(vbo_gl, GL_array_deleter {glDeleteBuffers});
    glBindBuffer(GL_ARRAY_BUFFER, render_data.vbo.get());
    glBufferData(
        GL_ARRAY_BUFFER,
        static_cast<GLsizei>(render_data.vertices.size() * sizeof(Vertex)),
        render_data.vertices.data(),
        GL_DYNAMIC_DRAW);

    GLuint ibo_gl {};
    glGenBuffers(1, &ibo_gl);
    render_data.ibo.reset(ibo_gl, GL_array_deleter {glDeleteBuffers});
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, render_data.ibo.get());
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                 static_cast<GLsizei>(render_data.indices.size() *
                                      sizeof(std::uint32_t)),
                 render_data.indices.data(),
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
}

void update_vertex_buffer(GLuint vao,
                          GLuint vbo,
                          const Render_data &render_data)
{
    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferSubData(
        GL_ARRAY_BUFFER,
        0,
        static_cast<GLsizei>(render_data.vertices.size() * sizeof(Vertex)),
        render_data.vertices.data());

    glBindVertexArray(0);
}

void create_render_geometry(Render_data &render_data)
{
    render_data.vertices.clear();
    render_data.indices.clear();

    render_data.line_indices_offset = render_data.indices.size();
    for (const auto &line : render_data.lines)
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
            static_cast<std::uint32_t>(render_data.vertices.size());
        render_data.vertices.push_back(
            {start_left,
             {-0.5f, 0.5f, -aspect_ratio - 0.5f, 0.0f},
             line.color});
        render_data.vertices.push_back(
            {start_right,
             {0.5f, 0.5f, -aspect_ratio - 0.5f, 0.0f},
             line.color});
        render_data.vertices.push_back(
            {end_right, {0.5f, -aspect_ratio - 0.5f, 0.5f, 0.0f}, line.color});
        render_data.vertices.push_back(
            {end_left, {-0.5f, -aspect_ratio - 0.5f, 0.5f, 0.0f}, line.color});
        render_data.indices.push_back(first_index + 0);
        render_data.indices.push_back(first_index + 1);
        render_data.indices.push_back(first_index + 2);
        render_data.indices.push_back(first_index + 0);
        render_data.indices.push_back(first_index + 2);
        render_data.indices.push_back(first_index + 3);
    }
    render_data.line_indices_size =
        render_data.indices.size() - render_data.line_indices_offset;

    render_data.circle_indices_offset = render_data.indices.size();
    for (const auto &circle : render_data.circles)
    {
        const auto half_side = circle.radius + 0.5f * circle.thickness;
        const auto bottom_left = circle.center + vec2 {-half_side, -half_side};
        const auto bottom_right = circle.center + vec2 {half_side, -half_side};
        const auto top_right = circle.center + vec2 {half_side, half_side};
        const auto top_left = circle.center + vec2 {-half_side, half_side};
        const auto rel_thickness = circle.thickness / half_side;

        const auto first_index =
            static_cast<std::uint32_t>(render_data.vertices.size());
        render_data.vertices.push_back(
            {bottom_left, {-1.0, -1.0f, rel_thickness, 0.0f}, circle.color});
        render_data.vertices.push_back(
            {bottom_right, {1.0f, -1.0f, rel_thickness, 0.0f}, circle.color});
        render_data.vertices.push_back(
            {top_right, {1.0f, 1.0f, rel_thickness, 0.0f}, circle.color});
        render_data.vertices.push_back(
            {top_left, {-1.0f, 1.0f, rel_thickness, 0.0f}, circle.color});
        render_data.indices.push_back(first_index + 0);
        render_data.indices.push_back(first_index + 1);
        render_data.indices.push_back(first_index + 2);
        render_data.indices.push_back(first_index + 0);
        render_data.indices.push_back(first_index + 2);
        render_data.indices.push_back(first_index + 3);
    }
    render_data.circle_indices_size =
        render_data.indices.size() - render_data.circle_indices_offset;
}

void init_application(Application &app)
{
    glfwSetErrorCallback(&glfw_error_callback);

    app.glfw_context.reset(glfwInit());
    CHECK(app.glfw_context.has_value());

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

    app.window.reset(
        glfwCreateWindow(1280, 720, "Truss Optimization", nullptr, nullptr));
    CHECK(app.window.has_value());

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
        centered_viewport(Application::world_size.x / Application::world_size.y,
                          app.framebuffer_width,
                          app.framebuffer_height);

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

    CHECK(IMGUI_CHECKVERSION());
    app.imgui_context.reset(ImGui::CreateContext());
    CHECK(app.imgui_context.has_value());

    auto &io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

    ImGui::StyleColorsDark();

    app.imgui_glfw_context.reset(
        ImGui_ImplGlfw_InitForOpenGL(app.window.get(), true));
    CHECK(app.imgui_glfw_context.has_value());

#ifdef __EMSCRIPTEN__
    ImGui_ImplGlfw_InstallEmscriptenCallbacks(app.window.get(), "#canvas");
#endif

    app.imgui_opengl_context.reset(ImGui_ImplOpenGL3_Init(glsl_version));
    CHECK(app.imgui_opengl_context.has_value());

    app.render_data.line_program = create_program(
        glsl_version, vertex_shader_code, line_fragment_shader_code);
    app.render_data.circle_program = create_program(
        glsl_version, vertex_shader_code, circle_fragment_shader_code);

    app.state = App_state::idle;

    app.problem = Problem::regular_grid;
    optimization_create_problem(app.optimization, app.problem);
}

void render_clear(Render_data &render_data)
{
    render_data.lines.clear();
    render_data.circles.clear();
}

void render_line(Render_data &render_data, const Line &line)
{
    render_data.lines.push_back(line);
}

void render_circle(Render_data &render_data, const Circle &circle)
{
    render_data.circles.push_back(circle);
}

void render(Render_data &render_data)
{
    create_render_geometry(render_data);
    create_vertex_and_index_buffers(render_data);

    glBindVertexArray(render_data.vao.get());

    glUseProgram(render_data.line_program.get());
    glDrawElements(GL_TRIANGLES,
                   static_cast<GLsizei>(render_data.line_indices_size),
                   GL_UNSIGNED_INT,
                   reinterpret_cast<void *>(render_data.line_indices_offset *
                                            sizeof(std::uint32_t)));

    glUseProgram(render_data.circle_program.get());
    glDrawElements(GL_TRIANGLES,
                   static_cast<GLsizei>(render_data.circle_indices_size),
                   GL_UNSIGNED_INT,
                   reinterpret_cast<void *>(render_data.circle_indices_offset *
                                            sizeof(std::uint32_t)));

    glBindVertexArray(0);
}

void make_ui(Application &app)
{
    if (ImGui::Begin("UI"))
    {
        ImGui::Text("%.3f ms/frame, %.2f fps",
                    1000.0 / static_cast<double>(ImGui::GetIO().Framerate),
                    static_cast<double>(ImGui::GetIO().Framerate));

        constexpr const char *items[] {"Regular grid", "Random Delaunay"};
        auto current_item = static_cast<int>(app.problem);
        ImGui::Combo("Problem", &current_item, items, std::size(items));
        app.problem = static_cast<Problem>(current_item);

        if (app.state == App_state::idle)
        {
            if (ImGui::Button("Start"))
            {
                app.state = App_state::running;
            }
        }
        else if (app.state == App_state::running)
        {
            if (ImGui::Button("Stop"))
            {
                app.state = App_state::idle;
            }
        }

        ImGui::SameLine();
        if (ImGui::Button("Step"))
        {
            app.state = App_state::running;
            app.should_idle = true;
        }

        ImGui::SameLine();
        if (ImGui::Button("Reset"))
        {
            app.step = 0;
            optimization_create_problem(app.optimization, app.problem);
        }
    }
    ImGui::End();
}

void main_loop_update(Application &app)
{
    glfwPollEvents();

    if (app.state == App_state::running)
    {
        optimization_step(app.optimization);
        ++app.step;

        if (app.should_idle)
        {
            app.state = App_state::idle;
            app.should_idle = false;
        }
    }

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    make_ui(app);

    ImGui::Render();

    render_clear(app.render_data);

    const auto [min_force, max_force] =
        std::ranges::minmax_element(app.optimization.axial_forces);

    for (std::size_t e {0}; e < app.optimization.elements.size(); ++e)
    {
        const auto [i, j] = app.optimization.elements[e];
        const auto activation = app.optimization.activations[e];

        vec3 color {1.0f, 1.0f, 1.0f};
        if (app.optimization.axial_forces.size() ==
            app.optimization.elements.size())
        {
            const auto force = app.optimization.axial_forces[e];
            const auto rel_force =
                force >= 0.0f ? force / *max_force : force / *min_force;
            const auto max_color = force >= 0.0f ? vec3 {0.25f, 0.25f, 1.0f}
                                                 : vec3 {1.0f, 0.25f, 0.25f};
            color = rel_force * max_color + 1.0f - rel_force;
        }

        render_line(app.render_data,
                    {.a = app.optimization.nodes[i],
                     .b = app.optimization.nodes[j],
                     .thickness = activation * 0.03f,
                     .color = color});
    }

    for (std::size_t i {0}; i < app.optimization.nodes.size(); ++i)
    {
        render_circle(app.render_data,
                      {.center = app.optimization.nodes[i],
                       .radius = 0.0f,
                       .thickness = 0.02f,
                       .color = {1.0f, 1.0f, 1.0f}});
    }

    // Draw gradient directions
    /*if (!app.optimization.gradients.empty())
    {
        // const auto max_gradient = *std::ranges::max_element(
        //     app.optimization.gradients,
        //     [](const vec2 &a, const vec2 &b) { return norm(a) < norm(b); });
        for (std::size_t i {0}; i < app.optimization.nodes.size(); ++i)
        {
            render_line(app.render_data,
                        {.a = app.optimization.nodes[i],
                         .b = app.optimization.nodes[i] -
                              0.07f * normalize(app.optimization.gradients[i]),
                         .thickness = 0.015f,
                         .color = {1.0f, 1.0f, 1.0f}});
        }
    }*/

    glViewport(app.viewport.x,
               app.viewport.y,
               app.viewport.width,
               app.viewport.height);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    render(app.render_data);

    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

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

    static Application app {};
    init_application(app);

    emscripten_set_main_loop_arg(
        [](void *arg) { main_loop_update(*static_cast<Application *>(arg)); },
        &app,
        0,
        true);

#else

    Application app {};
    init_application(app);

    while (!glfwWindowShouldClose(app.window.get()))
    {
        main_loop_update(app);
    }

#endif
}
