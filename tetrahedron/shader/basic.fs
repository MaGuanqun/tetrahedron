#version 460 core
out vec4 FragColor;

in VS_OUT{
	vec3 FragPos;
	vec3 Normal;
} fs_in;

struct Light{
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

uniform vec3 viewPos;
uniform vec3 lightPos;
uniform Light light;

uniform vec3 color;

uniform float transparence;


void main()
{	
    vec3 norm=normalize(fs_in.Normal);
	vec3 lightDir=normalize(lightPos-fs_in.FragPos);

    float diff=max(dot(norm,lightDir),0.0);
    vec3 diffuse=diff*light.diffuse;

	vec3 viewDir=normalize(viewPos-fs_in.FragPos);	
	vec3 reflectDir=reflect(-lightDir,norm);
	float spec=pow(max(dot(viewDir,reflectDir),0.0),32);
    vec3 specular=spec*light.specular;

    vec3 result=(light.ambient+diffuse+specular)*color;
    FragColor=vec4(result,transparence);

}
