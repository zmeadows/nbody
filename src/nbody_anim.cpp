#include "nbody_anim.h"

nbody_anim::nbody_anim(const size_t _width, const size_t _height)
    : width(_width), height(_height)
{
    solver = new nbody_solver(2,150,0.5);
}

nbody_anim::~nbody_anim(void ) 
{
    delete window;
    delete clock;
}

void nbody_anim::init_window(void) {
    window = new sf::Window(sf::VideoMode(width, height, 32), "SFML Window");
    clock = new sf::Clock;
    window->SetActive();
}

void nbody_anim::init_gl(void)
{
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho (0, width, height, 0, 0, 1);
    glDisable(GL_DEPTH_TEST);
    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(4.0);
    
    // Displacement trick for exact pixelization
    glTranslatef(0.375, 0.375, 0);
}

void nbody_anim::init(void )
{
    init_window();
    init_gl();
    loop();
}

void nbody_anim::loop(void )
{ 
  while (window->IsOpened()) {
    sf::Event event;
    while (window->GetEvent(event)) {
	// Close window : exit
	if (event.Type == sf::Event::Closed)
	    window->Close();

	// Escape key : exit
	if ((event.Type == sf::Event::KeyPressed) && (event.Key.Code == sf::Key::Escape))
	    window->Close();
    
	// Resize event : adjust viewport
	if (event.Type == sf::Event::Resized)
	    glViewport(0, 0, event.Size.Width, event.Size.Height);
    }
    
    
    if (clock->GetElapsedTime() > 0.02) {
	clock->Reset();
        glClearColor(.15, .15, .15, 0);
        glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_POINTS);
        size_t N = solver->get_num_bodies();
	for (size_t i = 0; i < N; ++i) {
            if (i == N-1) {
                glColor3f(0.9,0.6,0.2);
            } else { glColor3f(0.8,0.8,0.8); }
	    glVertex2f(5*solver->get_x_pos(i),4*solver->get_y_pos(i));
	}
	glEnd();
        window->Display();
    } else {
        solver->evolve_step();
    }
    

    
  }
}
