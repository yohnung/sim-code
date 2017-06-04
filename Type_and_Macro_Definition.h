//Type_and_Macro_Definition.h

// Logic value
#define True 1
#define False 0
// Type value
#define Incomplete 1      // Need to minus 1 from max grid number, 
                          // meaning every cycle starts at 0, ends at Grid_Num_?-1.
#define Complete 0        // Need not to minus 1.
// another Type value
#define Uniform 1         // Unifrom means that eta is setted to be a unifrom value
#define non_Unifrom 0    // non-Uniform means that eta can be abnormal value caused by, for example, fluctuation
// Order value
#define First 1
#define Second 2

#define Positive 1
#define Negative -1

typedef int Logic;        //Logic=True or False
typedef int Type;         // Type=Complete or Incomplete
typedef int Order;        // Order=First or Second
typedef int Symmetry_Type;// Symmetry_Type=Axial or Dot
