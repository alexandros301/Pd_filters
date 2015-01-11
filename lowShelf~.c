/****************************************************************
 * Low shelf filter external, combination of Mike Moser-Booth's *
 * [filtercoeff.mmb~]and [biquad.mmb~] abstractions             *
 * Written by Alexandros Drymonitis                             *
 ****************************************************************/

// Header files required by Pure Data
#include "m_pd.h"
#include "math.h"

// The class pointer
static t_class *lowShelf_class;

// The object structure
typedef struct _lowShelf {
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
} t_lowShelf;

// Function prototypes
void *lowShelf_new(void);
void lowShelf_dsp(t_lowShelf *x, t_signal **sp);
t_int *lowShelf_perform(t_int *w);

// The Pd class definition function
void lowShelf_tilde_setup(void)
{
	// Initialize the class
	lowShelf_class = class_new(gensym("lowShelf~"), (t_newmethod)lowShelf_new, 0, sizeof(t_lowShelf), 0, 0);

	// Specify signal input, with automatic float to signal conversion
	CLASS_MAINSIGNALIN(lowShelf_class, t_lowShelf, x_f);

	// Bind the DSP method, which is called when the DACs are turned on
	class_addmethod(lowShelf_class, (t_method)lowShelf_dsp, gensym("dsp"), A_CANT, 0);

	// Print authorship to Pd window
	post("lowShelf~: Filter external based on [filtercoeff.mmb~] and [biquad.mmb~] abstractions\nby Mike Moser-Booth, written by Alexandros Drymonitis");
}

// The new instance routine
void *lowShelf_new(void)
{
	// Basic object setup

	// Instantiate a new powSine~ object
	t_lowShelf *x = (t_lowShelf *) pd_new(lowShelf_class);

	// Create three additional singal inlets, the first one is on the house
	inlet_new(&x->obj, &x->obj.ob_pd, gensym("signal"), gensym("signal"));
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
t_int *lowShelf_perform(t_int *w)
{
	// Copy the object pointer
	t_lowShelf *x = (t_lowShelf *) (w[1]);

	// Copy signal vector pointers
	t_float *input = (t_float *) (w[2]);
	t_float *frequency = (t_float *) (w[3]);
	t_float *slope = (t_float *) (w[4]);
	t_float *gain = (t_float *) (w[5]);
	t_float *output = (t_float *) (w[6]);

	// Copy the signal vector size
	t_int n = w[7];

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
	float freq_local, slope_local, gain_local;
	float freq_cos, alpha;
	float in_sample, out_sample;

	// Perform the DSP loop
	while(n--){
		freq_local = *frequency++ * twopi / sr;
		slope_local = *slope++;
		gain_local = pow(10, *gain++/40);
		// clip slope's values
		if (slope_local < 0.001) slope_local = 0.001;
		else if (slope_local > 1) slope_local = 1;

		freq_cos = cos(freq_local);
		alpha = sin(freq_local) * sqrt((gain_local * gain_local + 1) * (1/slope_local - 1) + 2 * gain_local);

		// Coefficient values
		coeff_a0 = (gain_local + 1) + (gain_local - 1) * freq_cos + alpha;
		coeff_a1 = (-2) * ((gain_local - 1) + (gain_local + 1) * freq_cos);
		coeff_a2 = (gain_local + 1) + (gain_local - 1) * freq_cos - alpha;
		coeff_b0 = gain_local * ((gain_local + 1) - (gain_local - 1) * freq_cos + alpha);
		coeff_b1 = 2 * gain_local * ((gain_local - 1) - (gain_local + 1) * freq_cos);
		coeff_b2 = gain_local * ((gain_local + 1) - (gain_local - 1) * freq_cos - alpha);

		coeff_a1 = (coeff_a1 * (-1)) / coeff_a0;
		coeff_a2 = (coeff_a2 * (-1)) / coeff_a0;
		coeff_b0 /= coeff_a0;
		coeff_b1 /= coeff_a0;
		coeff_b2 /= coeff_a0;

		// This is copied from [biquad.mmb~]'s siggestion of [fexpr~] use
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
	return w + 8;
}

// The DSP method
void lowShelf_dsp(t_lowShelf *x, t_signal **sp)
{
	// Check if samplerate has changed
	if(x->x_sr != sp[0]->s_sr) x->x_sr = sp[0]->s_sr;

	/* Attach the object to the DSP chain, passing the DSP routine lowShelf_perform(),
	inlet and outlet pointers, and the signal vector size */
	dsp_add(lowShelf_perform, 7, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[0]->s_n);
}
