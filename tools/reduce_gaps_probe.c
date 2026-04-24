#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
    eGapAlignDel = 0,
    eGapAlignSub = 3,
    eGapAlignIns = 6
} EGapAlignOpType;

typedef struct {
    int size;
    EGapAlignOpType* op_type;
    int* num;
} GapEditScript;

static void s_UpdateEditScript(GapEditScript* esp, int pos, int bf, int af) {
    int op, qd, sd;

    if (bf > 0) {
        op = pos;
        qd = sd = bf;
        do {
            if (--op < 0) return;
            switch(esp->op_type[op]) {
            case eGapAlignSub:
                qd -= esp->num[op];
                sd -= esp->num[op];
                break;
            case eGapAlignIns:
                qd -= esp->num[op];
                break;
            case eGapAlignDel:
                sd -= esp->num[op];
            default:
                break;
            }
        } while (qd > 0 || sd > 0);

        esp->num[op] = -(qd > sd ? qd : sd);
        esp->op_type[op++] = eGapAlignSub;
        for (; op < pos-1; op++) esp->num[op] = 0;
        esp->num[pos] += bf;
        qd -= sd;
        esp->op_type[pos-1] = (qd>0) ? eGapAlignDel: eGapAlignIns;
        esp->num[pos-1] = (qd>0) ? qd : -qd;
    }

    if (af > 0) {
        op = pos;
        qd = sd = af;
        do {
            if (++op >= esp->size) return;
            switch(esp->op_type[op]) {
            case eGapAlignSub:
                qd -= esp->num[op];
                sd -= esp->num[op];
                break;
            case eGapAlignIns:
                qd -= esp->num[op];
                break;
            case eGapAlignDel:
                sd -= esp->num[op];
            default:
                break;
            }
        } while (qd > 0 || sd > 0);

        esp->num[op] = -(qd > sd ? qd : sd);
        esp->op_type[op--] = eGapAlignSub;
        for (; op > pos+1; op--) esp->num[op] = 0;
        esp->num[pos] += af;
        qd -= sd;
        esp->op_type[pos+1] = (qd>0) ? eGapAlignDel: eGapAlignIns;
        esp->num[pos+1] = (qd>0) ? qd : -qd;
    }
}

static void s_RebuildEditScript(GapEditScript* esp) {
    int i, j;
    for (i=0, j=-1; i<esp->size; i++) {
        if (esp->num[i] == 0) continue;
        if (j>=0 && esp->op_type[i] == esp->op_type[j]) {
            esp->num[j] += esp->num[i];
        } else if (j==-1 || esp->op_type[i] == eGapAlignSub
            || esp->op_type[j] == eGapAlignSub) {
            esp->op_type[++j] = esp->op_type[i];
            esp->num[j] = esp->num[i];
        } else {
            int d = esp->num[j] - esp->num[i];
            if (d > 0) {
                esp->num[j-1] += esp->num[i];
                esp->num[j] = d;
            } else if (d < 0) {
                if (j == 0 && i - j > 0) {
                    esp->op_type[j] = eGapAlignSub;
                    j++;
                }
                else {
                    esp->num[j-1] += esp->num[j];
                }
                esp->num[j] = -d;
                esp->op_type[j] = esp->op_type[i];
            } else {
                esp->num[j-1] += esp->num[j];
                --j;
            }
        }
    }
    esp->size = ++j;
}

static void s_ReduceGaps(GapEditScript* esp, const unsigned char *q, const unsigned char *s,
                         const unsigned char *qf, const unsigned char *sf) {
    int i, j, nm1, nm2, d;
    const unsigned char *q1, *s1;

    for (q1=q, s1=s, i=0; i<esp->size; i++) {
        if (esp->num[i] == 0) continue;
        if (esp->op_type[i] == eGapAlignSub) {
            if(esp->num[i] >= 12) {
                nm1 = 1;
                if (i > 0) {
                    while (q1-nm1>=q && (*(q1-nm1) == *(s1-nm1))) ++nm1;
                }
                q1 += esp->num[i];
                s1 += esp->num[i];
                nm2 = 0;
                if (i < esp->size -1) {
                    while ((q1+1<qf) && (s1+1<sf) && (*(q1++) == *(s1++))) ++nm2;
                }
                if (nm1>1 || nm2>0) s_UpdateEditScript(esp, i, nm1-1, nm2);
                q1--; s1--;
            } else {
                q1 += esp->num[i];
                s1 += esp->num[i];
            }
        } else if (esp->op_type[i] == eGapAlignIns) {
            q1 += esp->num[i];
        } else {
            s1 += esp->num[i];
        }
    }
    s_RebuildEditScript(esp);

    for (i=0; i<esp->size; i++) {
        if (esp->op_type[i] == eGapAlignSub) {
            q += esp->num[i];
            s += esp->num[i];
            continue;
        }
        if (i>1 && esp->op_type[i] != esp->op_type[i-2]
                && esp->num[i-2] > 0) {
            d = esp->num[i] + esp->num[i-1] + esp->num[i-2];
            if (d == 3) {
                (esp->num[i-2]) = 0;
                (esp->num[i-1]) = 2;
                (esp->num[i]) = 0;
                if (esp->op_type[i] == eGapAlignIns) {
                    ++q;
                } else {
                    ++s;
                }
            } else if (d < 12) {
                nm1 = 0;
                nm2 = 0;
                d = esp->num[i] < esp->num[i-2] ? esp->num[i] : esp->num[i-2];
                q -= esp->num[i-1];
                s -= esp->num[i-1];
                q1 = q;
                s1 = s;
                if (esp->op_type[i] == eGapAlignIns) {
                    s -= d;
                } else {
                    q -= d;
                }
                for (j=0; j<esp->num[i-1]; ++j, ++q1, ++s1, ++q, ++s) {
                    if (*q1 == *s1) nm1++;
                    if (*q == *s) nm2++;
                }
                for (j=0; j<d; ++j, ++q, ++s) {
                    if (*q == *s) nm2++;
                }
                if (nm2 >= nm1 - d) {
                    (esp->num[i-2]) -= d;
                    (esp->num[i-1]) += d;
                    (esp->num[i]) -= d;
                } else {
                    q = q1;
                    s = s1;
                }
            }
        }
        if (esp->op_type[i] == eGapAlignIns) {
            q += esp->num[i];
        } else {
            s += esp->num[i];
        }
    }
    s_RebuildEditScript(esp);
}

static int hexval(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

static unsigned char* parse_hex(const char* hex, int* len) {
    int n = (int)strlen(hex) / 2;
    unsigned char* out = (unsigned char*)malloc((size_t)n);
    for (int i = 0; i < n; i++) {
        out[i] = (unsigned char)((hexval(hex[2*i]) << 4) | hexval(hex[2*i+1]));
    }
    *len = n;
    return out;
}

static EGapAlignOpType parse_op(const char* op) {
    if (strcmp(op, "Ins") == 0) return eGapAlignIns;
    if (strcmp(op, "Del") == 0) return eGapAlignDel;
    return eGapAlignSub;
}

int main(void) {
    char qhex[8192], shex[8192], op[16];
    int nops;
    if (scanf("%8191s %8191s %d", qhex, shex, &nops) != 3) return 2;
    GapEditScript esp;
    esp.size = nops;
    esp.op_type = (EGapAlignOpType*)malloc(sizeof(EGapAlignOpType) * (size_t)nops);
    esp.num = (int*)malloc(sizeof(int) * (size_t)nops);
    for (int i = 0; i < nops; i++) {
        if (scanf("%15s %d", op, &esp.num[i]) != 2) return 3;
        esp.op_type[i] = parse_op(op);
    }
    int qlen, slen;
    unsigned char* q = parse_hex(qhex, &qlen);
    unsigned char* s = parse_hex(shex, &slen);
    s_ReduceGaps(&esp, q, s, q + qlen, s + slen);
    for (int i = 0; i < esp.size; i++) {
        const char* name = esp.op_type[i] == eGapAlignIns ? "Ins" :
                           esp.op_type[i] == eGapAlignDel ? "Del" : "Sub";
        printf("%s %d%s", name, esp.num[i], i + 1 == esp.size ? "\n" : " ");
    }
    free(q);
    free(s);
    free(esp.op_type);
    free(esp.num);
    return 0;
}
