#include <math.h>
#include <stdio.h>
#include "cnx.h"
#include "cnxbasic.h"
#include "cnxpa.h"
#include "cnxrand.h"
#include "cnxadapt.h"

#include "iir.h"

/*
typedef struct
{
	int t;
	CNX_Connector *inbound_ref;
	CNX_Connector *inbound_sig;
} GainPrinter;

// processor only safe for queue bound connectors
int ProcessGainPrinter(void *gainprinter, int count)
{
	int n;
	int *t = &((GainPrinter*)gainprinter)->t;
	CNX_Connector *inbound_ref = ((GainPrinter*)gainprinter)->inbound_ref;
	CNX_Connector *inbound_sig = ((GainPrinter*)gainprinter)->inbound_sig;
	int bufferlength = inbound_ref->bufferlength;
	for(n = 0; n < count; n++)
	{
		CNX_Buffer *in_ref = CNX_TakeQueue(inbound_ref);
		CNX_Buffer *in_sig = CNX_TakeQueue(inbound_sig);
		float *ref = in_ref->data;
		float *sig = in_sig->data;
		float refp = 0;
		float sigp = 0;
		for(int i = 0; i < bufferlength; i++)
		{
			refp += ref[i] * ref[i];
			sigp += sig[i] * sig[i];
		}
		if(((*t)++)%16==0) printf("Gain: %f %f %fdB\n", refp, sigp, (10*log10(sigp/refp)));
		CNX_PoolBuffer(inbound_ref, in_ref);
		CNX_PoolBuffer(inbound_sig, in_sig);
	}
	return n;
}
*/

int main(int argc, char **argv)
{
	IIR_FilterSpec spec;
	spec.design = IIR_BUTTERWORTH;
	spec.Rp = 1;
	spec.As = 10;
	/*
	spec.response = IIR_LOWPASS;
	spec.freq.lp.dw = 0.1;
	spec.freq.lp.wc = 2;
	*/
	/*
	spec.response = IIR_HIGHPASS;
	spec.freq.hp.dw = 0.1;
	spec.freq.hp.wc = 2;
	*/

	spec.response = IIR_BANDPASS;
	spec.freq.bp.dw = 0.1;
	spec.freq.bp.wl = 1.5;
	spec.freq.bp.wu = 2.0;

	IIR_DesignFilter(&spec);
	/*
	CNXPA_Initialize();
	CNX_Connector noise_splitter1;
	CNX_Connector splitter1_amplifier;
	CNX_Connector amplifier_stream;
	CNX_Connector stream_splitter2;
	CNX_Connector splitter1_adaptive;
	CNX_Connector splitter2_adaptive;
	CNX_Connector splitter2_printer;
	CNX_Connector adaptive_printer;
	CNXRAND_Generator noise;
	CNXBASIC_Splitter splitter1;
	CNXBASIC_Amplifier amplifier;
	CNXPA_Stream stream;
	CNXBASIC_Splitter splitter2;
	CNXADAPT_Adaptive adaptive;
	GainPrinter printer;
	int length = 128;
	CNX_ConstructConnector(&noise_splitter1, sizeof(CNX_Scalar), length, 4, CNX_POOLBOUND | CNX_QUEUEBOUND);
	CNX_ConstructConnector(&splitter1_amplifier, sizeof(CNX_Scalar), length, 0, CNX_POOLBOUND | CNX_QUEUEBOUND);
	CNX_ConstructConnector(&splitter1_adaptive, sizeof(CNX_Scalar), length, 0, CNX_POOLBOUND | CNX_QUEUEBOUND);
	CNX_ConstructConnector(&amplifier_stream, sizeof(CNX_Scalar), length, 4, CNX_POOLBOUND | CNX_QUEUEBOUND);
	CNX_ConstructConnector(&stream_splitter2, sizeof(CNX_Scalar), length, 4, CNX_QUEUEBOUND);
	CNX_ConstructConnector(&splitter2_adaptive, sizeof(CNX_Scalar), length, 0, CNX_QUEUEBOUND);
	CNX_ConstructConnector(&splitter2_printer, sizeof(CNX_Scalar), length, 0, CNX_QUEUEBOUND);
	CNX_ConstructConnector(&adaptive_printer, sizeof(CNX_Scalar), length, 0, CNX_QUEUEBOUND);
	CNXRAND_ConstructGenerator(&noise, &noise_splitter1, 42);
	splitter1.inbound = &noise_splitter1;
	splitter1.outbounds = malloc(3*sizeof(CNX_Connector*));
	splitter1.outbounds[0] = &splitter1_amplifier;
	splitter1.outbounds[1] = &splitter1_adaptive;
	splitter1.outbounds[2] = NULL;
	amplifier.gain = 0.1;
	amplifier.inbound = &splitter1_amplifier;
	amplifier.outbound = &amplifier_stream;
	CNXPA_ConstructStream(&stream, &amplifier_stream, &stream_splitter2, 8000, length);
	splitter2.inbound = &stream_splitter2;
	splitter2.outbounds = malloc(3*sizeof(CNX_Connector*));
	splitter2.outbounds[0] = &splitter2_adaptive;
	splitter2.outbounds[1] = &splitter2_printer;
	splitter2.outbounds[2] = NULL;
	CNXADAPT_ConstructAdaptive(&adaptive, length, 512, 0.01, 1.0, &splitter1_adaptive, &splitter2_adaptive, NULL, &adaptive_printer, NULL);
	printer.inbound_ref = &splitter2_printer;
	printer.inbound_sig = &adaptive_printer;
	CNX_Instance *noiseinst = CNX_StartInstance(CNXRAND_Generate, &noise);
	CNX_Instance *splitter1inst = CNX_StartInstance(CNXBASIC_ProcessRoundRobin, &splitter1);
	CNX_Instance *amplifierinst = CNX_StartInstance(CNXBASIC_ProcessAmplifier, &amplifier);
	CNX_Instance *splitter2inst = CNX_StartInstance(CNXBASIC_ProcessRoundRobin, &splitter2);
	CNX_Instance *adaptiveinst = CNX_StartInstance(CNXADAPT_ProcessAdaptive, &adaptive);
	CNX_Instance *printerinst = CNX_StartInstance(ProcessGainPrinter, &printer);
	//CNX_Instance *consumer1inst = CNX_StartInstance(CNXBASIC_ProcessConsumer, &robin1_adaptive);
	//CNX_Instance *consumer2inst = CNX_StartInstance(CNXBASIC_ProcessConsumer, &robin2_adaptive);
	//CNX_Instance *consumer3inst = CNX_StartInstance(CNXBASIC_ProcessConsumer, &adaptive_printer);
	//CNX_Instance *consumer4inst = CNX_StartInstance(CNXBASIC_ProcessConsumer, &robin2_printer);
	CNXPA_StartStream(&stream);
	getchar();
	printf("stopping\n");
	CNXPA_StopStream(&stream);
	printf("audio stopped\n");
	CNXPA_Terminate();
	printf("pa terminated\n");
	return 0;
	*/
}
