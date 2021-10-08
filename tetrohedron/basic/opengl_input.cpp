#include"opengl_input.h"
#include "../external/glfw3.h"
#define CONTAINS(list,element) ((list).find(element) != (list).end())
#define REMOVE(list,element) ((list).erase((list).find(element)))

void MouseInput::beginFrame() {
	memcpy(last_screen_pos, screen_pos, 16);

	prev_left_press = left_press;
	prev_right_press = right_press;

	mouse_callback = false;
	scroll_callback = false;
}

void KeyboardInput::beginFrame() {
	previous = current;
}

void KeyboardInput::callback(int key, int action) {
	if (action == GLFW_PRESS) {
		current.insert(key);
	}
	else if (action == GLFW_RELEASE) {
		if (CONTAINS(current, key)) REMOVE(current, key);
	}
}

bool KeyboardInput::keyIsPressed(int key) const { return (CONTAINS(current, key)); }
bool KeyboardInput::keyWasPressedThisFrame(int key) const { return (CONTAINS(current, key) && !CONTAINS(previous, key)); }
bool KeyboardInput::keyWasReleasedThisFrame(int key) const { return (!CONTAINS(current, key) && CONTAINS(previous, key)); }



void Input::beginFrame() {

	keyboard.beginFrame();
	mouse.beginFrame();
}

void Input::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (!guiCaptureKeyboard) {
		keyboard.callback(key, action);
	}
}


void Input::mouseCallback(GLFWwindow* window, double xpos, double ypos) {
	mouse.screen_pos[0] = xpos;
	mouse.screen_pos[1] = ypos;

	mouse.mouse_callback = true;

	if (mouse.right_press || mouse.left_press) {
		mouse.angle[0] = M_PI * (xpos- mouse.last_screen_pos[0]) / SCR_WIDTH;
		mouse.angle[1] = M_PI * (mouse.last_screen_pos[1] - ypos) / SCR_HEIGHT;
		mouse.move_direction[0]= xpos - mouse.last_screen_pos[0];
		mouse.move_direction[1] = mouse.last_screen_pos[1] - ypos;
	}
}

void Input::mouseButtonCallback(GLFWwindow* window, int button, int action, int mod)
{
	if (!guiCaptureMouse) {
		int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
		if (state == GLFW_PRESS) {
			mouse.left_press = true;
		}
		else if (state == GLFW_RELEASE) {
			mouse.left_press = false;
		}

		int right_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
		if (right_state == GLFW_PRESS) {
			mouse.right_press = true;
		}
		else if (right_state == GLFW_RELEASE) {
			mouse.right_press = false;
		}
	}
}

void Input::scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
	if (!guiCaptureMouse) {
		mouse.scroll_callback = true;
		mouse.scroll = yoffset;
	}
}