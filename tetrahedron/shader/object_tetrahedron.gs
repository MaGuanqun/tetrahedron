#version 460 core
layout (triangles) in;
layout (triangle_strip, max_vertices=3) out;

in vec3 FragPos1[];
out vec3 Normal;
out vec3 FragPos;

void main()
{
    vec3 a = (FragPos1[1] - FragPos1[0]).xyz;
    vec3 b = (FragPos1[2] - FragPos1[0]).xyz;
    Normal = normalize(cross(a, b));
    for(int i = 0; i < gl_in.length(); ++i)
    {

        gl_Position = gl_in[i].gl_Position;
        FragPos=FragPos1[i];
        EmitVertex();
    }
    EndPrimitive();
  
} 