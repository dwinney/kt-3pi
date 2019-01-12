
#ifndef _AMPLITUDE_
#define _AMPLITUDE_

class Decay
{
public:

const int qn_J, qn_C, qn_P, qn_I, qn_H;   // quantum numbers of decaying particle

void setJPC(int j, int p, int c);
void setIsospin(int i);
void setHelicity(int lambda);
};

class Isobar : Amplitude
{
public:
};

#endif
