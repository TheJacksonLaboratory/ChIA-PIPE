#ifndef CPU_ADAPTER_H
#define CPU_ADAPTER_H

// |5' adapter| = 33
#define G_NEXTERA_5PAdapter "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
// |3' adapter| = 34
#define G_NEXTERA_3PAdapter "GTCTCGTGGGCTCGGAGATGTGTATAAGAGAGAG"

#define ADAPTER_NONE   0
#define ADAPTER_5P     1
#define ADAPTER_3P     2
#define ADAPTER_5Prc   4
#define ADAPTER_3Prc   8

static char adapterTypeBuffer[11];
static const char *adapterTypeToStr(int32_t t) {
	if (ADAPTER_NONE==t) return ".";
	else if (ADAPTER_5P==t) return "5P";
	else if (ADAPTER_3P==t) return "3P";
	else if (ADAPTER_5Prc==t) return "rc5P";
	else if (ADAPTER_3Prc==t) return "rc3P";
	else {
		//return "?";
		sprintf(adapterTypeBuffer, "%#x", t);
		return adapterTypeBuffer;
	}
}

static int trim_read_3prime(int trim_qual, int minReadLen, bseq1_t *seqs)
{
	int s = 0, l, max = 0, max_l = seqs->l_seq;
	if (trim_qual < 1) return max_l;
	for (l = seqs->l_seq - 1; l >= minReadLen; --l) {
		s += trim_qual - (seqs->qual[l] - 33);
		if (s < 0) break;
		if (s > max) max = s, max_l = l;
	}
	return max_l;
}

#endif
