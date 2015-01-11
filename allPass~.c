/***************************************************************
 * All pass filter external, combination of Mike Moser-Booth's *
 * [filtercoeff.mmb~] abstraction and Pd's [biquad~]           *
 * Written by Alexandros Drymonitis                            *
 ***************************************************************/

// Header files required by Pure Data
#include "m_pd.h"
#include "math.h"

// The class pointer
static t_class *allPass_class;

// The object structure
typedef struct _allPass {
	// The Pd object
	t_object obj;
	// Convert floats to signals
	t_float x_f;
	// Filter coefficient variables
	float x_coeff_a0, x_coeff_a1, x_coeff_a2;
	float x_coeff_b0, x_coeff_b1, x_coeff_b2;
	// Previous samples variables
	float x_last_in, x_prev_in, x_last_out, x_prev_out;
	// Rest of variables
	float x_twopi;
	float x_sr;
} t_allPass;

// Function prototypes
void *allPass_new(void);
void allPass_dsp(t_allPass *x, t_signal **sp);
t_int *allPass_perform(t_int *w);

// The Pd class definition function
void allPass_tilde_setup(void)
{
	// Initialize the class
	allPass_class = class_new(gensym("allPass~"), (t_newmethod)allPass_new, 0, sizeof(t_allPass), 0, 0);

	// Specify signal input, with automatic float to signal conversion
	CLASS_MAINSIGNALIN(allPass_class, t_allPass, x_f);

	// Bind the DSP method, which is called when the DACs are turned on
	class_addmethod(allPass_class, (t_method)allPass_dsp, gensym("dsp"), A_CANT, 0);

	// Print authorship to Pd window
	post("allPass~: Filter external based on [filtercoeff.mmb~] and [biquad.mmb~] abstractions\nby Mike Moser-Booth,  written by Alexandros Drymonitis");
}

// The new instance routine
void *allPass_new(void)
{
	// Basic object setup

	// Instantiate a new powSine~ object
	t_allPass *x = (t_allPass *) pd_new(allPass_class);

	// Create two additional singal inlets, the first one is on the house
	inlet_new(&x->obj, &x->obj.ob_pd, gensym("signal"), gensym("signal"));
	inlet_new(&x->obj, &x->obj.ob_pd, gensym("signal"), gensym("signal"));

	// Create one signal outlet
	outlet_new(&x->obj, gensym("signal"));

	// Initialize coefficient variables to zero
	x->x_coeff_a0 = 0;
	x->x_coeff_a1 = 0;
	x->x_coeff_a2 = 0;
	x->x_coeff_b0 = 0;
	x->x_coeff_b1 = 0;
	x->x_coeff_b2 = 0;

	// Initialize previous samples variables to zero
	x->x_last_in = 0;
	x->x_prev_in = 0;
	x->x_last_out = 0;
	x->x_prev_out = 0;

	// define twoPi
	x->x_twopi = 8.0 * atan(1.0);

	// get system's sampling rate
	x->x_sr = sys_getsr();

	// Return a pointer to the new object
	return x;
}

// The perform routine
t_int *allPass_perform(t_int *w)
{
	// Copy the object pointer
	t_allPass *x = (t_allPass *) (w[1]);

	// Copy signal vector pointers
	t_float *input = (t_float *) (w[2]);
	t_float *frequency = (t_float *) (w[3]);
	t_float *bandwidth = (t_float *) (w[4]);
	t_float *output = (t_float *) (w[5]);

	// Copy the signal vector size
	t_int n = w[6];

	// Dereference components from the object structure
	float coeff_a0 = x->x_coeff_a0;
	float coeff_a1 = x->x_coeff_a1;
	float coeff_a2 = x->x_coeff_a2;
	float coeff_b0 = x->x_coeff_b0;
	float coeff_b1 = x->x_coeff_b1;
	float coeff_b2 = x->x_coeff_b2;
	float last_in = x->x_last_in;
	float prev_in = x->x_prev_in;
	float last_out = x->x_last_out;
	float prev_out = x->x_prev_out;
	float twopi = x->x_twopi;
	float sr = x->x_sr;
	// Local variables
	float freq_local, bandwidth_local;
	float freq_sin, alpha;
	float in_sample, out_sample;

	// Perform the DSP loop
	while(n--){
		freq_local = *frequency++ * twopi / sr;
		bandwidth_local = *bandwidth++;
		bandwidth_local = (bandwidth_local > 0.01 ? bandwidth_local : 0.01);
		freq_sin = sin(freq_local);
		alpha = freq_sin * sinh(log(2) / 2 * bandwidth_local * freq_local / freq_sin);

		// Coefficient values
		coeff_a0 = alpha + 1;
		coeff_a1 = (-2) * cos(freq_local);
		coeff_a2 = 1 - alpha;
		coeff_b0 = 1 - alpha;
		coeff_b1 = (-2) * cos(freq_local);
		coeff_b2 = alpha + 1;

		coeff_a1 = (coeff_a1 * (-1)) / coeff_a0;
		coeff_a2 = (coeff_a2 * (-1)) / coeff_a0;
		coeff_b0 /= coeff_a0;
		coeff_b1 /= coeff_a0;
		coeff_b2 /= coeff_a0;

		// This is copied from [biquad.mmb~]'s suggestion of [fexpr~] use
		in_sample = *input++;
		out_sample = (coeff_b0 * in_sample) + (coeff_b1 * last_in) + (coeff_b2 * prev_in) + (coeff_a1 * last_out) + (coeff_a2 * prev_out);
		*output++ = out_sample;
		prev_in = last_in;
		last_in = in_sample;
		prev_out = last_out;
		last_out = out_sample;
		}

	// Update object's coefficients
	x->x_coeff_a0 = coeff_a0;
	x->x_coeff_a1 = coeff_a1;
	x->x_coeff_a2 = coeff_a2;
	x->x_coeff_b0 = coeff_b0;
	x->x_coeff_b1 = coeff_b1;
	x->x_coeff_b2 = coeff_b2;

	// Update object's previous samples
	x->x_last_in = last_in;
	x->x_prev_in = prev_in;
	x->x_last_out = last_out;
	x->x_prev_out = prev_out;

	// Return the next address in the DSP chain
	return w + 7;
}

// The DSP method
void allPass_dsp(t_allPass *x, t_signal **sp)
{
	// Check if samplerate has changed
	if(x->x_sr != sp[0]->s_sr) x->x_sr = sp[0]->s_sr;

	/* Attach the object to the DSP chain, passing the DSP routine allPass_perform(),
	inlet and outlet pointers, and the signal vector size */
	dsp_add(allPass_perform, 6, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
}
