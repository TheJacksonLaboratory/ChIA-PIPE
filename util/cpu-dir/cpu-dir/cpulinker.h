#ifndef CPU_LINKER_H
#define CPU_LINKER_H

#define G_CP_HalfLinkerA "GTTGGATAAGATATCGCGG"
#define G_CP_HalfLinkerB "GTTGGAATGTATATCGCGG"

#define LINKER_NONE         0x0000

#define LINKER_PAIRED       0x0010 // paired tags
#define LINKER_TIED         0x0020 // tied linker
#define LINKER_HALF         0x0040 // a half linker
#define LINKER_FULL         0x0080 // a full linker

#define LINKER_RC_A_PRIV    0x0001 // right most linkerA = rc(linkerA)
#define LINKER_A_PRIV       0x0002 // left most linkerA
#define LINKER_RC_B_PRIV    0x0004 // right most linkerB = rc(linkerB)
#define LINKER_B_PRIV       0x0008 // left most linkerB
#define LINKER_MASK         0x000F // to retrieve all combination of linker

#define LINKER_RC_A   (LINKER_RC_A_PRIV|LINKER_HALF)
#define LINKER_A      (LINKER_A_PRIV|LINKER_HALF)
#define LINKER_RC_B   (LINKER_RC_B_PRIV|LINKER_HALF)
#define LINKER_B      (LINKER_B_PRIV|LINKER_HALF)

#define LINKER_AA     (LINKER_A_PRIV|LINKER_RC_A_PRIV|LINKER_FULL)
#define LINKER_BB     (LINKER_B_PRIV|LINKER_RC_B_PRIV|LINKER_FULL)
#define LINKER_AB     (LINKER_A_PRIV|LINKER_RC_B_PRIV|LINKER_FULL)
#define LINKER_BA     (LINKER_B_PRIV|LINKER_RC_A_PRIV|LINKER_FULL)

#define LINKER_AA_PRIV (LINKER_A_PRIV|LINKER_RC_A_PRIV)
#define LINKER_BB_PRIV (LINKER_B_PRIV|LINKER_RC_B_PRIV)
#define LINKER_AB_PRIV (LINKER_A_PRIV|LINKER_RC_B_PRIV)
#define LINKER_BA_PRIV (LINKER_B_PRIV|LINKER_RC_A_PRIV)

static inline int isFullLinker(int l) {
	return (l&LINKER_FULL);
}

static inline int isHalfLinker(int l) {
	return (l&LINKER_HALF);
}

static inline int isNonChimericFullLinker(int l) {
	return (l==LINKER_AA || l==LINKER_BB);
}

static inline int isTiedLinker(int l) {
	return (l&LINKER_TIED);
}

static inline int isPairedLinker(int l) {
	return (l&LINKER_PAIRED);
}

// TODO: not thread-safe
//4,294,967,296
static char linkerTypeBuffer[11];
static const char *linkerTypeToStr(int16_t t) {
	if (LINKER_NONE==t) return ".";
	else if (LINKER_AA==(LINKER_AA&t)) return "Aa";
	else if (LINKER_BB==(LINKER_BB&t)) return "Bb";
	else if (LINKER_AB==(LINKER_AB&t)) return "Ab";
	else if (LINKER_BA==(LINKER_BA&t)) return "Ba";
	else if (LINKER_A==(LINKER_A&t)) return "A";
	else if (LINKER_B==(LINKER_B&t)) return "B";
	else if (LINKER_RC_A==(LINKER_RC_A&t)) return "a";
	else if (LINKER_RC_B==(LINKER_RC_B&t)) return "b";
	else {
		//return "?";
		sprintf(linkerTypeBuffer, "%#x", t);
		return linkerTypeBuffer;
	}
}

static char shortLinkerTypeBuffer[11];
static const char *shortLinkerTypeToStr(int16_t t) {
	int16_t maskt = t & LINKER_MASK;
	if (LINKER_NONE==maskt) return ".";
	else if (LINKER_AA_PRIV==(LINKER_AA_PRIV&maskt)) return "Aa";
	else if (LINKER_BB_PRIV==(LINKER_BB_PRIV&maskt)) return "Bb";
	else if (LINKER_AB_PRIV==(LINKER_AB_PRIV&maskt)) return "Ab";
	else if (LINKER_BA_PRIV==(LINKER_BA_PRIV&maskt)) return "Ba";
	else if (LINKER_A_PRIV==(LINKER_A_PRIV&maskt)) return "A";
	else if (LINKER_B_PRIV==(LINKER_B_PRIV&maskt)) return "B";
	else if (LINKER_RC_A_PRIV==(LINKER_RC_A_PRIV&maskt)) return "a";
	else if (LINKER_RC_B_PRIV==(LINKER_RC_B_PRIV&maskt)) return "b";
	else {
		//return "?";
		sprintf(shortLinkerTypeBuffer, "%#x", t);
		return shortLinkerTypeBuffer;
	}
}

static char linkerClassBuffer[11];
static const char *linkerClassToStr(int16_t t) {
	if (LINKER_NONE==t) return "none";
	else if (LINKER_PAIRED==(LINKER_PAIRED&t)) return "paired";
	else if (LINKER_TIED==(LINKER_TIED&t)) return "tied";
	else if (LINKER_HALF==(LINKER_HALF&t)) return "half";
	else if (LINKER_FULL==(LINKER_FULL&t)) return "full";
	else {
		//return "?";
		sprintf(linkerClassBuffer, "%#x", t);
		return linkerClassBuffer;
	}
}

#endif
