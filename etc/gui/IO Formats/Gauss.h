#pragma once
#include "IOSpinSystem.h"

public ref class Gauss : SpinachGUI::IOSpinSystem
{
public:
    Gauss(System::String^ filename, Model^% model);
    virtual bool LoadFile() override;
};
