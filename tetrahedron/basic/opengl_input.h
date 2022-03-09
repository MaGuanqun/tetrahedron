#pragma once
#include <unordered_set>
#include"global.h"
#include"../external/glad.h"
#include"../external/glfw3.h"
// Handles mouse button input
struct MouseInput {
	double last_screen_pos[2] = { 0,0 };
	double screen_pos[2] = { 0.0, 0.0 };

	float angle[2] = { 0.0,0.0 };
	double move_direction[2] = { 0.0,0.0 };

	bool left_press = false;
	bool prev_left_press = false;
	bool right_press = false;
	bool prev_right_press = false;

	double scroll = 0.0;

	bool mouse_callback = false;
	bool scroll_callback = false;

	// Must be called once per frame. Typically handled by Input class.
	void beginFrame();

	bool leftButtonIsPressed() { return left_press; }
	bool leftButtonWasPressedThisFrame() { return left_press && !prev_left_press; }
	bool leftButtonWasReleasedThisFrame() { return !left_press && prev_left_press; }
	bool leftButtonWasPressedPreviousAndThisFrame() { return left_press && prev_left_press; }

	bool rightButtonIsPressed() { return right_press; }
	bool rightButtonWasPressedThisFrame() { return right_press && !prev_right_press; }
	bool rightButtonWasReleasedThisFrame() { return !right_press && prev_right_press; }
};


struct KeyboardInput {
	KeyboardInput() = default;

	KeyboardInput(const KeyboardInput&) = delete;
	KeyboardInput& operator=(const KeyboardInput&) = delete;

	// Must be called once per frame, before callback
	void beginFrame();
	void callback(int key, int action);

	bool keyIsPressed(int key) const;
	bool keyWasPressedThisFrame(int key) const;
	bool keyWasReleasedThisFrame(int key) const;
private:
	// Set of the keys that were being pressed this frame
	std::unordered_set<int> current;
	// Set of the keys that were being pressed the previous frame
	std::unordered_set<int> previous;
};

struct Input {
	MouseInput mouse;
	KeyboardInput keyboard;

	Input() = default;
	Input(const Input&) = delete;
	Input& operator=(const Input&) = delete;
	void beginFrame();
	void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
	void mouseCallback(GLFWwindow* window, double xpos, double ypos);
	void mouseButtonCallback(GLFWwindow* window, int button, int action, int mod);
	void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);

	bool guiCaptureMouse = false;
	bool guiCaptureKeyboard = false;
};