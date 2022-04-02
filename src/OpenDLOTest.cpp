// DLOTest.cpp
// Copyright (c) 2022, Renier Cloete
// This program is released under the BSD license. See the file LICENSE.txt for details.

#include "Domain.h"
#include "Enums.h"
#include "MosekDLOSolver.h"
#include "CoinDLOSolver.h"

#include <glad/gl.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include "linmath.h"

#include <stdlib.h>
#include <stdio.h>

static const char* vertex_shader_text =
"#version 110\n"
"uniform mat4 MVP;\n"
"attribute vec3 vCol;\n"
"attribute vec2 vPos;\n"
"varying vec3 color;\n"
"void main()\n"
"{\n"
"    gl_Position = MVP * vec4(vPos, 0.0, 1.0);\n"
"    color = vCol;\n"
"}\n";

static const char* fragment_shader_text =
"#version 110\n"
"varying vec3 color;\n"
"void main()\n"
"{\n"
"    gl_FragColor = vec4(color, 1.0);\n"
"}\n";

static void error_callback( int error, const char* description )
{
	fprintf( stderr, "Error: %s\n", description );
}

static void key_callback( GLFWwindow* window, int key, int scancode, int action, int mods )
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose( window, GLFW_TRUE );
}

int main()
{
	//Calculate yield lines & load factor

	CDomain Domain;

	Domain.AddBoundaryPoint( { 0.0, 0.0 }, eEdgeType::FIXED );
	Domain.AddBoundaryPoint( { 1.0, 0.0 }, eEdgeType::FIXED );
	Domain.AddBoundaryPoint( { 1.0, 1.0 }, eEdgeType::FIXED );
	Domain.AddBoundaryPoint( { 0.0, 1.0 }, eEdgeType::FIXED );
	
	Domain.Discretize( 0.25 );
	Domain.BuildEdges();
	Domain.SetLoads( 1, 0 );
	Domain.SetYieldMoments( 1, 1, 1, 1 );

	std::vector<double> YieldEdges;
	std::vector<CPoint2D> DomainEdges;
	
	DomainEdges = Domain.GetBoundaryPoints();

	CCoinDLOSolver Solver;
	double Lambda = Solver.Solve( &Domain );

	//Output

	std::cout << "Lambda: " << Lambda;

	YieldEdges = Solver.GetEdgeData();

	struct sVertex
	{
		float x, y;
		float r, g, b;
	};
	std::vector<sVertex> vertices;

	auto addVtx = [&vertices]( float x, float y, float mx, float my, float m, float r, float g, float b)
	{
		vertices.push_back( { x / m - mx / 2, y / m - my / 2, r, g, b } );
	};

	float max_x = -1;
	float max_y = -1;
	float max = -1;

	if (DomainEdges.size())
	{
		CPoint2D p1, p2;

		for (size_t i = 0; i < DomainEdges.size(); ++i)
		{
			p1 = DomainEdges[i];

			if (p1.x > max_x) max_x = (float)p1.x;
			if (p1.y > max_y) max_y = (float)p1.y;
		}

		max = std::max( max_x, max_y );
		max_x = max_x / max;
		max_y = max_y / max;

		for (size_t i = 0; i + 1 < DomainEdges.size(); ++i)
		{
			p1 = DomainEdges[i];
			p2 = DomainEdges[i + 1];

			addVtx( (float)p1.x, (float)p1.y, max_x, max_y, max, 0.0f, 1.0f, 0.0f );
			addVtx( (float)p2.x, (float)p2.y, max_x, max_y, max, 0.0f, 1.0f, 0.0f );
		}

		p1 = DomainEdges[DomainEdges.size()-1];
		p2 = DomainEdges[0];

		addVtx( (float)p1.x, (float)p1.y, max_x, max_y, max, 0.0f, 1.0f, 0.0f );
		addVtx( (float)p2.x, (float)p2.y, max_x, max_y, max, 0.0f, 1.0f, 0.0f );
	}

	for (size_t i = 0; i < YieldEdges.size() / 8; ++i)
	{
		size_t i1 = 8 * i;
		double pm = YieldEdges[i1 + 3];

		float r = 0.0f, g = 0.0f, b = 0.0f;

		if (abs( pm ))
		{
			if (YieldEdges[i1] > 0)
				r = 1.0f;
			else
			{
				g = 1.0f;
				b = 1.0f;
			}

			addVtx( (float)YieldEdges[i1 + 4], (float)YieldEdges[i1 + 5], max_x, max_y, max, r, g, b );
			addVtx( (float)YieldEdges[i1 + 6], (float)YieldEdges[i1 + 7], max_x, max_y, max, r, g, b );
		}
	}

	GLFWwindow* window;
	GLuint vertex_buffer, vertex_shader, fragment_shader, program;
	GLint mvp_location, vpos_location, vcol_location;

	glfwSetErrorCallback( error_callback );

	if (!glfwInit())
		exit( EXIT_FAILURE );

	glfwWindowHint( GLFW_CONTEXT_VERSION_MAJOR, 2 );
	glfwWindowHint( GLFW_CONTEXT_VERSION_MINOR, 0 );

	window = glfwCreateWindow( 640, 480, "OpenDLO example", NULL, NULL );
	if (!window)
	{
		glfwTerminate();
		exit( EXIT_FAILURE );
	}

	glfwSetKeyCallback( window, key_callback );

	glfwMakeContextCurrent( window );
	gladLoadGL( glfwGetProcAddress );
	glfwSwapInterval( 1 );

	glGenBuffers( 1, &vertex_buffer );
	glBindBuffer( GL_ARRAY_BUFFER, vertex_buffer );
	glBufferData( GL_ARRAY_BUFFER, vertices.size()*sizeof( sVertex ), &vertices[0], GL_STATIC_DRAW);

	vertex_shader = glCreateShader( GL_VERTEX_SHADER );
	glShaderSource( vertex_shader, 1, &vertex_shader_text, NULL );
	glCompileShader( vertex_shader );

	fragment_shader = glCreateShader( GL_FRAGMENT_SHADER );
	glShaderSource( fragment_shader, 1, &fragment_shader_text, NULL );
	glCompileShader( fragment_shader );

	program = glCreateProgram();
	glAttachShader( program, vertex_shader );
	glAttachShader( program, fragment_shader );
	glLinkProgram( program );

	mvp_location = glGetUniformLocation( program, "MVP" );
	vpos_location = glGetAttribLocation( program, "vPos" );
	vcol_location = glGetAttribLocation( program, "vCol" );

	glEnableVertexAttribArray( vpos_location );
	glVertexAttribPointer( vpos_location, 2, GL_FLOAT, GL_FALSE,
						   sizeof( vertices[0] ), (void*)0 );
	glEnableVertexAttribArray( vcol_location );
	glVertexAttribPointer( vcol_location, 3, GL_FLOAT, GL_FALSE,
						   sizeof( vertices[0] ), (void*)(sizeof( float ) * 2) );

	while (!glfwWindowShouldClose( window ))
	{
		float ratio;
		int width, height;
		mat4x4 m, p, mvp;

		glfwGetFramebufferSize( window, &width, &height );
		ratio = width / (float)height;

		glViewport( 0, 0, width, height );
		glClear( GL_COLOR_BUFFER_BIT );

		mat4x4_identity( m );
		mat4x4_ortho( p, -ratio, ratio, -1.f, 1.f, 1.f, -1.f );
		mat4x4_mul( mvp, p, m );

		glUseProgram( program );
		glUniformMatrix4fv( mvp_location, 1, GL_FALSE, (const GLfloat*)mvp );
		glDrawArrays( GL_LINES, 0, (GLsizei)vertices.size() );

		glfwSwapBuffers( window );
		glfwPollEvents();
	}

	glfwDestroyWindow( window );

	glfwTerminate();
	exit( EXIT_SUCCESS );

}
