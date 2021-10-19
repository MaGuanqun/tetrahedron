#version 460 core                                                              

out vec3 FragColor;
uniform int picking_base_ID;
uniform vec4 xy_range;
void main()                                                                         
{                  
   int objectID = gl_PrimitiveID+picking_base_ID;
   int r,g,b;
   int resi;
   b=objectID/65536;
   resi=objectID%65536;
   g=resi/256;
   r=resi%256;
   if(gl_FrontFacing){
	   FragColor = vec3(float(r) / 255.0f, 
	  float(g) / 255.0f, 
	   float(b) / 255.0f);      
  }
  else{
	FragColor=vec3(1.0,1.0,1.0);
  }
}