#pragma once
#include "IOSpinSystem.h"

using namespace SpinachGUI;

public ref class Orca : IOSpinSystem
{
public:
    Orca(String^ filename, Model^% model);
    virtual bool LoadFile() override;

private:
    ::Matrix3x3^ ReadMatrix();
    void AddTensor(Atom^ a, Atom^ b, InteractionKind kind, ::Matrix3x3^ matrix, Unit^ unit);
};
