import regex
import sys
import time

from argparse import ArgumentParser
from collections import defaultdict


REPEATS = {}
COMPL = dict(zip("ACGT", "TGCA"))
def rc(s):
    return "".join(list(map(lambda x: x in COMPL and COMPL[x] or x, s))[::-1])


def warn(x):
    sys.stderr.write(str(x)+"\n")

IUPAC_WILDCARDS = {
    "W" : "AT",
    "S" : "CG",
    "M" : "AC",
    "K" : "GT",
    "R" : "AG",
    "Y" : "CT",
    "B" : "CGT",
    "D" : "AGT",
    "H" : "ACT",
    "V" : "ACG",
    "N" : "ACGT",
}

def re_disamb(seq):
    out = ""
    wildcards = frozenset(IUPAC_WILDCARDS.keys())
    for c in seq:
        if c in wildcards:
            out += "[" + "|".join(IUPAC_WILDCARDS[c]) + "]"
        else:
            out += c
    return out

def get_parts(repeat):
    seq = repeat.replace(" ","")
    has_o = "(" in seq
    has_c = ")" in seq
    out_seq = seq.replace("(","").replace(")","")
    out_before, out_after = None, None
    if has_o and has_c:
        out_before = seq[seq.find("(")+1:].replace("(","").replace(")","")
        out_after = seq[:seq.find(")")].replace("(","").replace(")","")
    return (out_seq, out_before, out_after)

def load_repeats(filename):
    out = defaultdict(dict)
    warn("\nloading repeats from " + filename)
    with open(filename, "r") as file:
        for line in file:
            cp = line.find("#")
            if cp != -1:
                line = line[0:cp]
            line = line.strip()
            if not line:
                continue
            if line.count("\t") != 2:
                warn("malformed line: \n" + line + "\nskipping...")
                continue
            sys, repeat_raw, mm_pre = line.split("\t")
            mm, ins, dl, _ = (mm_pre+",0,0,0,0").split(",",3)
            repeat, opn, cls = get_parts(repeat_raw)
            out[sys].update({ repeat : { "MM" : int(mm), "INS" : int(ins), "DEL" : int(dl) } })
            if opn and cls:
                out[sys][repeat].update({"OPEN" : opn, "CLOSE" : cls})
    return out or None

def fill_parts(frm, all_mm, id2primer = None):
    for sys in frm.keys():
        for seq in frm[sys].keys():
            mm = frm[sys][seq]["MM"]
            ins = frm[sys][seq]["INS"]
            dl = frm[sys][seq]["DEL"]
            if "OPEN" in frm[sys][seq] and "CLOSE" in frm[sys][seq]:
                before_seq = frm[sys][seq]["OPEN"]
                after_seq = frm[sys][seq]["CLOSE"]
            else:
                coord = len(seq)/2
                if coord < 18:
                    coord = 18
                before_seq = seq[-coord:]
                after_seq = seq[:coord]

            for s in [seq, before_seq, after_seq]:
                recalc = False
                if not s in all_mm:
                    all_mm[s] = { "MM" : mm , "INS" : ins, "DEL" : dl, "SYS" : set(), "TYPE" : set(), "ID" : str(len(all_mm)), "DA_SET" : set() }
                    recalc = True
                if mm > all_mm[s]["MM"]:
                    all_mm[s]["MM"] = mm
                    recalc = True
                if recalc:
                    all_mm[s]["SYS"].add(sys)
                    # overlapped here or in find?
                    regex_s = re_disamb(s)
                    regex_str = "(?:%s){s<=%d,i<=%d,d<=%d}" % (regex_s, mm, ins, dl)
                    #regex_str = "(?:%s){s<=%d}" % (s, mm)
                    all_mm[s]["RE_FUZZY"] = regex.compile(regex_str, overlapped = True)
                if s == seq:
                    all_mm[s]["TYPE"].update([sys+"(",sys+")"])
                elif s == before_seq:
                    all_mm[s]["TYPE"].update([sys+"("])
                elif s == after_seq:
                    all_mm[s]["TYPE"].update([sys+")"])

    if id2primer is not None  and type(id2primer) == dict:
        for seq in all_mm:
            seq_id = all_mm[seq]["ID"]
            id2primer.update({seq_id: seq})

    warn("\nusing following repeats sequences:")
    for seq in all_mm:
        dump = "seq: " + seq + " id: " + all_mm[seq]["ID"] + " types: " + "|".join(all_mm[seq]["TYPE"]) + " mismatches: " + str(all_mm[seq]["MM"])
        warn(dump)

def mm_dist(s1, s2, max_dist = -1):
    l1 = len(s1)
    l2 = len(s2)
    mm = l2 - l1
    if mm < 0:
        mm = -mm
    check_dist = max_dist >= 0
    for i in range(min(l1, l2)):
        if s1[i] != s2[i]:
            mm += 1
        if check_dist and mm > max_dist:
            return -1
    return mm

def fuzzy_find_seq(where, what, mm, ins, dl, what_re_fuzzy, pat_type = None, pat_id = None):
    '''returns [(mm, from, to, substr)]'''
    if mm < 0:
        mm = 0
    out = []
    positions = [(i.span()[0], i.span()[1], i[0], i.fuzzy_counts[0], i.fuzzy_counts[1], i.fuzzy_counts[2]) for i in what_re_fuzzy.finditer(where, overlapped=True)]
    positions = filter(lambda x: x[3] <= mm and x[4] <= ins and x[5] <= dl, positions)
    for frm, to, match, dist_mm, dist_ins, dist_del in positions:
        out.append((dist_mm + dist_ins + dist_del, frm, to, match, pat_type, pat_id))
    if not out:
        return None
    return out

def prune_matches(matches, sys):
    # [(1, 207L, 225L, 'GTTGCAAGGGATTGGGCC'), (1, 207L, 243L, 'GTTGCAAGGGATTGGGCCCCGTAAGGGGATTGCGAC')
    # join halves
    joined = []
    prev = None
    #print "\nsys", sys, "\n"
    # TODO: FIX MM NUMBER
    for t in sorted(matches, key = lambda x: (x[1],x[2])):
        #print t
        if prev and t[1] != prev[1]:
            if t[2] != prev[2]:
                joined.append(prev)
                prev = t
        else:
            prev = t
    if prev:
        joined.append(prev)
    # filter for open close pairs, dublicate if necessary
    sys_o, sys_c = sys+"(", sys+")"
    sys_o_c_set = set([sys_o, sys_c])
    res = []
    prev = None
    for t in joined:
       #print t
       if not prev:
           prev = t
           continue
       t_type = t[4].intersection(sys_o_c_set)
       prev_type = prev[4]
       # append only open-close pairs
       if sys_o in prev_type and sys_c in t_type:
          res.append(prev)
          res.append(t)
       prev = t
    if len(res) < 2:
        return None
    return res or None


def get_system_coords(seq, all_mm):
    out = defaultdict(list)
    for pat in all_mm:
        mm = all_mm[pat]["MM"]
        ins = all_mm[pat]["INS"]
        dl = all_mm[pat]["DEL"]
        pat_sys = all_mm[pat]["SYS"]
        pat_type = all_mm[pat]["TYPE"]
        pat_re_fuzzy = all_mm[pat]["RE_FUZZY"]
        pat_id = all_mm[pat]["ID"]
        found = fuzzy_find_seq(seq, pat, mm, ins, dl, pat_re_fuzzy, pat_type, pat_id)
        if found:
            for sys in pat_sys:
                out[sys] += found
    res = {}
    for sys in out:
        pruned = prune_matches(out[sys], sys)
        if pruned:
            res[sys] = pruned
    return res


def str_time(start = 0):
    delta = int(time.time() - start)
    return ":".join(map(lambda i: "%02d" % i, [ delta/3600, (delta/60) % 60, delta % 60]))

def diff_str(orig, diff):
    return "".join(map(lambda o,d: o and d and o!=d and d.lower() or d or "", orig.upper(), diff.upper()))
    #return diff

def get_prm_parts(tpl, id2prm):
    mm, f, t, seq, st, seq_id = tpl
    seq_orig = seq_id in id2prm and id2prm[seq_id] or ""
    if len(st) == 2:
        st = ")("
    else:
        st = list(st)[0][-1]
    return (f, t, seq_orig, mm, diff_str(seq_orig, seq), st)


def repeat_counts(arr):
    size = len(arr)
    for cont, cont_size in enumerate(arr):
        for i in range(cont_size):
            yield size, cont, cont_size, i

def dump_read(place, read_id, res, seq, id2prm, both_dirs = False, tag = "", qual = "", dump_linked_ids = False):
    if not res:
        return
    systems_cnt = len(res)
    systems = ",".join(sorted(res.keys()))
    overlap = "F"
    if systems_cnt != 1:
        prev = None
        for p in sorted([ (crds[0][1], crds[-1][2]) for crds in res.values() ]):
            if prev:
                if prev[1] > p[0]:
                    overlap = "T"
                    break
            prev = p

    both_dirs = both_dirs and "T" or "F"

    for sys in res:
        coords = res[sys]
        if len(coords) % 2 != 0:
            warn(" ".join(map(str(place, read_id, tag, sys, "coords not paired... skipping last not paired"))))
        total = len(coords) / 2

        # deal with adjacent primers and zero-sized spacers
        linkage_size = [0]
        # assume that len % 2 == 0 and adding fake end
        for start, end, new_start in zip(coords[0::2], coords[1::2], coords[2::2] + [coords[-1]]):
            open_from, open_to, open_orig, mm, open_observed, open_type = get_prm_parts(start, id2prm)
            close_from, close_to, close_orig, mm, close_observed, close_type = get_prm_parts(end, id2prm)
            spacer_len = close_from - open_to
            # nospacer case
            if spacer_len <= 0:
                if linkage_size[-1] != 0:
                    linkage_size.append(0)
                linkage_size[-1] += 1
                linkage_size.append(0)
                continue
            linkage_size[-1] += 1
            # different opening sequence
            if end != new_start:
                linkage_size.append(0)
        # removing trailing 0 for empty or ending with nospacer region
        if linkage_size[-1] == 0:
             linkage_size.pop()

        for i, pair_cont in enumerate(zip(coords[0::2], coords[1::2], repeat_counts(linkage_size))):
             open_from, open_to, open_orig, mm, open_observed, open_type = get_prm_parts(pair_cont[0], id2prm)
             close_from, close_to, close_orig, mm, close_observed, close_type = get_prm_parts(pair_cont[1], id2prm)
             spacer_len = close_from - open_to
             spacer = seq[open_to:close_from]
             cont_dump = []
             if dump_linked_ids:
                 cont_dump = list(pair_cont[2])
             qual_dump = []
             if qual:
                 qual_dump = [
                     qual[open_from:open_to],
                     qual[close_from:close_to],
                     qual[open_to:close_from]
                 ]
             print ("\t".join(map(str, [
                place, read_id, tag, both_dirs, systems_cnt, systems, overlap, sys,
                open_from, open_to, open_orig, mm, open_observed, open_type,
                close_from, close_to, close_orig, mm, close_observed, close_type,
                total, i, spacer_len, spacer
             ]+ qual_dump + cont_dump)))


def main():
    # debug
    # keep_cryptic
    # min_primer_len [18]
    # keep_non_overlapping_systems

    # trim N at the beginning and end(?) (different length though)
    # keep_edge_ns
    # i-e.1
    # remove those with overlapping systems
    # fwd and rev?
    # discard if there's other sys primers inbetween
    parser = ArgumentParser()
    parser.add_argument("--primers", required=True, help="tab-separated 3 cols primers file")
    parser.add_argument("--each", required=False, help="each line of how many to process, format '0/4', for parallel runs")
    parser.add_argument("--time-each", required=False, help="report times for every TIME_EACH processed(and skipped) lines", type=int, default=50000)
    parser.add_argument("--time-first", required=False, help="report times after TIME_FIRST initially processed / skipped lines", type=int, default=1000)
    parser.add_argument("--dump-linked-ids", required=False, help="appends spacer ids withing linked regions (# linked regions, id of region, size of region, spacer id)", action="store_true")
    args = parser.parse_args()

    part, of_parts = None, None
    try:
        part, of_parts = map(int, args.each.strip().split("/", 1))
    except:
        pass
    if of_parts:
        warn("\nid %d of %d, zero-started" %(part, of_parts))

    repeats = load_repeats(args.primers)
    if not repeats:
        warn("no repeats loaded, exiting")
        exit(1)

    all_mm = {}
    id2primers = {}
    fill_parts(repeats, all_mm, id2primers)

    warn("\ntrimming...")

    start_time = time.time()
    for cnt, line in enumerate(sys.stdin, start = 1):
        if cnt % args.time_each == 0 or cnt == args.time_first:
            warn("\nprocessed " + str(cnt) + " reads in " + str_time(start_time) + " ...")
        if of_parts and cnt % of_parts != part:
            continue

        rec = list(map(lambda x: x.strip(), line.strip().split("\t")))
        rec += [""]
        place, read_id, seq, qual = rec[0:4]

        rc_seq = rc(seq)
        rc_qual = qual[::-1]

        res = get_system_coords(seq, all_mm)
        res_rc = get_system_coords(rc_seq, all_mm)

        both_dirs = res and res_rc
        dump_read(place, read_id, res, seq, id2primers, both_dirs, "+", qual, args.dump_linked_ids)
        dump_read(place, read_id, res_rc, rc_seq, id2primers, both_dirs, "-", qual, args.dump_linked_ids)

    warn("\nprocessed " + str(cnt) + " reads in " + str_time(start_time) + " ...")
    warn("\ndone...")


if __name__ == "__main__":
    main()

