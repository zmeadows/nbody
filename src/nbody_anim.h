#ifndef NBODY_ANIM_H
#define NBODY_ANIM_H

#include "nbody.h"
#include <SFML/Window.hpp>

class nbody_anim {
private:
    const size_t width;
    const size_t height;
    
    nbody_solver *solver;
    
    sf::Window *window;
    sf::Clock *clock;
    
    void init_window(void);
    void init_gl(void);
    void end_window(void);
    void end_gl(void);
    
    void loop(void);
public:
    nbody_anim(
	const size_t _width,
	const size_t _height
    );
    ~nbody_anim(void);
    
    void init(void);
};

#endif
