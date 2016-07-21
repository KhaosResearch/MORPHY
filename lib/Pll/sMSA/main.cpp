#include <iostream>

#include "Menus.h"

using namespace std;

int main(int argc, char* argv[]) {
    // Check the number of parameters
    if (argc < 2) {
        Menus* menu = new Menus();
        menu->ShowMainMenu();
    }else{
        Menus* menu = new Menus();
        menu->InitTripleInitializer(1, false);
        menu->ExecuteResampling();
    }

    return 0;
}
