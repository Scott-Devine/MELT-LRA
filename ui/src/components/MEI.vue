<script setup>

  import { reactive } from 'vue'
  import { color } from 'd3-color'
  import { interpolateYlOrRd } from 'd3-scale-chromatic'
  import { format } from 'd3-format'
  import { scaleLinear, scaleLog } from 'd3-scale';
  import { axisTop, axisBottom, axisLeft, axisRight } from 'd3-axis';
  
  const props = defineProps({
    meiStr: { type: String, required: true },
    label: { type: String, required: false }
  })

  const ME_LENGTHS = {
    'ALU': 281,
    'SVA': 1316,
    'SVA_A': 1387,
    'SVA_F': 1375,
    'LINE1': 6019
  }

  const csv_headers = ['chrom', 'pos', 'strand', 'ME', '%ME', '%id', '%id_ng', '%cov', 'insertion_seq', 'left_flank_seq', 'right_flank_seq', 'TSD_seq', 'polyX_coords', 'ME_coords', 'insertion_coords', 'match_string']

  function parse_mei(mei_str) {
    const f = mei_str.split(',')
    let nh = csv_headers.length
    let mei = {};
    for (let i = 0;i < nh; ++i) {
      mei[csv_headers[i]] = f[i];
    }
    mei['pos'] = mei['pos'] * 1.0
    return mei;
  }
//  const mei_str = "chr22,16414896,+,LINE1,3.3%,34.2%,90.4%,97.1%,ATTTAGAGACTTGTCCAAGATGTATATTAGTTCCTTGACTCTCTGTCCTAACGTAGTCCAGAGACCTGAGCTGTCTGGACATTCTGTAGCATGTTACATTTGTCCATTTT,TTCACATTCATGGGAAGGACAACAGCATGC,ATTTACCAGTTAAAATAGACTGGGTAAGAA,ATTTA,111-110,678-876,5-107,|.||||^||.|vvvv||||^^^^^^^^|.|.|||.||||||^|||vvvvvvvvvvv.||||||vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv|||.|||v|.||||vvvvvv|||||vvvv||||||vvvvvvvvvvvvvvvvvvvvvvvvvvv||||||vvvvvvvvvvvvvvvv||||||^|||^^^^||^^^^^||||"
//  const mei_str = "chr22,11860400,+,ALU,100.0%,97.5%,97.9%,100.0%,AAAATGAATGTATACGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAAACCATCCCGGCTAAAACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAAATTAGCCGGGCGTAGTGGCGGGCGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCCCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAA,TTCCCAGTTACATGGGATCTATAGTTCTAC,AAAATGAATGTATACATAATCCATGTAAAT,AAAATGAATGTATAC,298-312,1-281,16-297,||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.|||||||.||||||.||||||||||||||||||||||||||||||^||||||||||||||||||.||||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||||||||||||||"
//  const mei_str = "chr22,17224400,-,ALU,100.0%,94.1%,96.8%,100.0%,AACAAGTGCTAATAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTGAGACGGAGTCTCGCTCTGTCGCCCAGGCCGGACTGCGGACTGCAGTGGCGCAATCTCGGCTCACTGCAAGCTCCGCTTCCCGGGTTCACGCCATTCTCCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCGGCCACCGCGCCCGGCTAATTTTTTTGTATTTTTAGTAGAGACGGGGTTTCACCTTGTTAGCCAGGATGGTCTCGATCTCCTGACCTCATGATCCACCCGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACCGCGCCCGGCC,AAAAAAAACTTGGAAGATGGAAGGTACATA,AACAAGTGCTAATAATTTAGAAGACAAAAA,AACAAGTGCTAATAATTT,19-44,1-281,45-333,||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||||.|||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||||^|||||||||||||||||.||||||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||||||||||||||||||||||.||||||||^^^^^^^|||.|||.|||||||||||||||||||||||||||||"
//  const mei = parse_mei(mei_str)
 
  const mei = parse_mei(props.meiStr)
  
  // SVG coordinate system
  const width = 1200
  const height = 250
  const margins = { left: 80, right: 75, top: 40, bottom: 20}
  let rx1 = margins.left
  let rx2 = width - margins.right
 
  const left_flank_bp = mei.left_flank_seq.length
  const right_flank_bp = mei.right_flank_seq.length
  const flank_bp = left_flank_bp + right_flank_bp

  // genomic coordinate system: insertion + flanking seq
  const ins_len = mei.insertion_seq.length
  const ref_xscale = scaleLinear().domain([0, ins_len + flank_bp]).range([rx1, rx2])
  
  // insertion coordinate system
  const ins_rx1 = ref_xscale(left_flank_bp)
  const ins_rx2 = ref_xscale(left_flank_bp + ins_len)
  const ins_xscale = scaleLinear().domain([0, ins_len]).range([ins_rx1, ins_rx2]);
  const ins_xaxis = axisTop(ins_xscale)
  const[px_x1, px_x2] = mei.polyX_coords.split("-").map(x => x * 1.0)
  const[ins_x1, ins_x2] = mei.insertion_coords.split("-").map(x => x * 1.0)

  // genomic flanking regions
  const lflank_xscale = scaleLinear().domain([mei.pos-left_flank_bp, mei.pos]).range([rx1, ins_rx1])
  const lflank_xaxis = axisTop(lflank_xscale).tickValues([mei.pos]).tickFormat(bp => mei.chrom + ":" + bp).tickSizeOuter(0)
  const rflank_xscale = scaleLinear().domain([mei.pos + ins_len, mei.pos + ins_len + right_flank_bp]).range([ins_rx2, rx2])
  const rflank_xaxis = axisTop(rflank_xscale).tickValues([]).tickFormat(bp => "").tickSizeOuter(0)

  // reference ME coordinate system
  const me_len = ME_LENGTHS[mei.ME]
 
  // srhink ME to fit
  let me_rx1 = ins_rx1
  let me_rx2 = ins_rx2
  
  // but use same scale if ME is smaller
  if (me_len <= ins_len + flank_bp) {
    const me_diff = (ins_len + flank_bp) - me_len
    me_rx1 = ref_xscale(me_diff/2)
    me_rx2 = ref_xscale(me_diff/2 + me_len)
  } 

  const me_domain = mei.strand == '+' ? [0, me_len] : [me_len, 0]
  const me_xscale = scaleLinear().domain(me_domain).range([me_rx1, me_rx2]);
  const me_xaxis = axisBottom(me_xscale)
  const[me_x1, me_x2] = mei.ME_coords.split("-").map(x => x * 1.0)

  // convert match string to alignment spans
  let spans = [];
  // offsets from left side of the alignment - will add ins_x1, me_x1 to these
  let ins_o1 = 0
  let ins_o2 = 0
  let me_o1 = 0
  let me_o2 = 0
  let n_id_bp = 0
  const msl = mei.match_string.length;

  function add_span() {
    if (((ins_o2 - ins_o1) != 0) && ((me_o2 - me_o1) != 0)) {
      // handle reverse strand matches
      let ins_c = null 
      if (mei.strand == '+') {
        ins_c = [ins_x1 + ins_o1 - 1, ins_x1 + ins_o2 - 1]
      } else {
        ins_c = [ins_x2 - ins_o1 - 1, ins_x2 - ins_o2 - 1]
      }
      const span = {
        'ins': ins_c, 
        'me': [me_x1 + me_o1 - 1, me_x1 + me_o2 - 1],
        'pct_id': (n_id_bp / (ins_o2 - ins_o1)) * 100.0
      }
      spans.push(span)
      n_id_bp = 0
    }
  }

  // set color based on percent identity
  const color_fn = interpolateYlOrRd;
  let col = null;

  function pctid_color(pctid, alpha) {
    const clr = color(color_fn(pctid / 100.0))
    clr.opacity = alpha
    return clr.formatRgb()
  }

  // match string contains only ^ (gap in insertion), v (gap in ME), |, and .ß
  for(let i = 0;i < msl; ++i) {
    if (mei.match_string[i] == '^') {
      add_span()
      ins_o2 += 1
      ins_o1 = ins_o2
      me_o1 = me_o2
    } else if (mei.match_string[i] == 'v') {
      add_span()
      me_o2 += 1
      ins_o1 = ins_o2
      me_o1 = me_o2
    } else {
      if (mei.match_string[i] == '|') n_id_bp += 1
      // start or extend span
      ins_o2 += 1
      me_o2 += 1
    }
  }
  add_span()

  const feat_y_offset = 6
  const lbl_y_offset = 22
  const ins_axis_y = margins.top + 50
  const me_axis_y = height - margins.bottom - 50

  // convert spans to polygons
  const match_y1 = ins_axis_y + feat_y_offset
  const match_y2 = me_axis_y - feat_y_offset

  spans.forEach(s => {
    s.points = []
    s.points.push([ins_xscale(s.ins[0]), match_y1])
    s.points.push([ins_xscale(s.ins[1]), match_y1])
    s.points.push([me_xscale(s.me[1]), match_y2])
    s.points.push([me_xscale(s.me[0]), match_y2])
    s.points_str = s.points.join(" ")
  })

  const me_ref_str = ": " + mei['%ME'] + " coverage at " + mei['%id_ng'] + " identity"

  // color key
  const cb_height = 20
  const cb_width = 15
  const cb_y = height - margins.bottom - (cb_height * 5)
  const color_key_blocks = [
    { 'id': 95, 'color': pctid_color(95, 1), 'y': cb_y, height: cb_height },
    { 'id': 90, 'color': pctid_color(90, 1), 'y': cb_y + cb_height, 'height': cb_height },
    { 'id': 80, 'color': pctid_color(80, 1), 'y': cb_y + cb_height * 2, 'height': cb_height },
    { 'id': 70, 'color': pctid_color(70, 1), 'y': cb_y + cb_height * 3, 'height': cb_height },
    { 'id': 60, 'color': pctid_color(60, 1), 'y': cb_y + cb_height * 4, 'height': cb_height },
  ]

  const state = reactive({
    // SVG
    width: width,
    height: height,
    margins: margins,

    // props
    label: props.label,
    // MEI
    mei: mei, 
    ins_len: ins_len,
    tsd_len: mei.TSD_seq.length,
    // flanking regions
    rflank_x: mei.pos + ins_len,
    lflank_xscale: lflank_xscale,
    lflank_xaxis: lflank_xaxis,
    rflank_xscale: rflank_xscale,
    rflank_xaxis: rflank_xaxis,
    // insertion and associated features
    ins_xscale: ins_xscale, 
    ins_xaxis: ins_xaxis,
    polyx_x1: px_x1,
    polyx_x2: px_x2,
    ins_x1: ins_x1,
    ins_x2: ins_x2,
    // reference ME
    me_len: me_len,
    me_xscale: me_xscale,
    me_xaxis: me_xaxis,
    me_x1: me_x1,
    me_x2: me_x2,
    // y offsets
    // insertion feature on ref genome axis (top)
    ins_label_y: margins.top - 10,
    ins_y: margins.top + 3,
    // axes for left and right genomic flanking regions (top)
    flank_axis_y: margins.top,
    flank_feat_y: margins.top + feat_y_offset,
    flank_feat_lbl_y: margins.top + lbl_y_offset,
    // insertion axis and features (below genomic axis)
    ins_axis_y: ins_axis_y,
    ins_feat_y: ins_axis_y + feat_y_offset,
    ins_feat_lbl_y: ins_axis_y + lbl_y_offset,
    // mobile element reference 
    me_axis_y: me_axis_y,
    me_lbl_y: me_axis_y + 50,
    // alignment spans
    spans: spans,
    me_ref_str: me_ref_str,
    // color key
    color_key_blocks: color_key_blocks
  })
  
</script>

<template>
  <div class="greetings pb-4">
    <h3 class="green">{{state.label}} [{{state.mei.ME}}] {{ state.mei.chrom }}:{{ state.mei.pos }}:{{ state.mei.strand }}:{{state.ins_len}} bp</h3>

    <svg
      ref="mei_svg"
      key="mei-svg-1"
      :width="state.width"
      :height="state.height"
      class="pa-0 ma-0"
      style="border: 1px solid white"
      >

      <!-- ad-hoc color key -->
      <rect v-for="ck in color_key_blocks" :x="10" :y="ck.y" :width="20" :height="ck.height" :fill="ck.color"></rect>
      <text v-for="ck in color_key_blocks" :x="35" :y="ck.y + ck.height -5" fill="#d0d0d0">{{ck.id + "%"}}</text>

      <!-- insertion feature on ref genome -->
      <line :x1="state.ins_xscale(0)" :x2="state.ins_xscale(state.ins_len)" :y1="state.ins_y" :y2="state.ins_y" stroke-width="3" stroke="#ffffff" />
      <text class="ins_label" :x="state.ins_xscale(state.ins_len/2)" :y="state.ins_label_y" fill="#ffffff">{{ state.ins_len + " bp insertion [" + state.mei.strand + " strand]" }}</text>
      <!-- genomic sequence flanking regions -->
      <g v-axis="state.lflank_xaxis" class="xaxis" :transform="`translate(0,${state.flank_axis_y})`"></g>
      <g v-axis="state.rflank_xaxis" class="xaxis" :transform="`translate(0,${state.flank_axis_y})`"></g>

      <!-- TSD in right flanking sequence-->
      <line v-if="state.tsd_len > 0" :x1="state.rflank_xscale(state.rflank_x)" :x2="state.rflank_xscale(state.rflank_x + state.tsd_len)" :y1="state.flank_feat_y" :y2="state.flank_feat_y" stroke-width="4" stroke="#a0ffa0"/>
      <text v-if="state.tsd_len > 0" :x="state.rflank_xscale(state.rflank_x)" :y="state.flank_feat_lbl_y" fill="#ffffff">TSD</text>

      <!-- inserted sequence axis -->
      <g v-axis="state.ins_xaxis" class="xaxis" :transform="`translate(0,${state.ins_axis_y})`">
      </g>

      <!-- reference ME axis -->
      <g v-axis="state.me_xaxis" class="xaxis" :transform="`translate(0,${state.me_axis_y})`">
      </g>
      <text v-if="state.mei.strand == '+'" class="me_label" :x="state.me_xscale(0)" :y="state.me_lbl_y" fill="#ffffff" text-anchor="start">{{ state.me_len + " bp " + state.mei.ME + " reference " + me_ref_str + " >>>"}}</text>
      <text v-else class="me_label" :x="state.me_xscale(0)" :y="state.me_lbl_y" fill="#ffffff" text-anchor="end">{{ "<<< " + state.me_len + " bp " + state.mei.ME + " reference " + me_ref_str}}</text>

      <!-- match/alignment spans -->
      <polygon v-for="s in state.spans" :points="s.points_str" :fill="pctid_color(s.pct_id, 0.5)" :stroke="pctid_color(s.pct_id, 1.0)" stroke-width="2" />

      <!-- TSD, polyX in ME coords -->
      <line v-if="state.tsd_len > 0" :x1="state.ins_xscale(0)" :x2="state.ins_xscale(state.tsd_len)" :y1="state.ins_feat_y" :y2="state.ins_feat_y" stroke-width="4" stroke="#a0ffa0"/>
      <text v-if="state.tsd_len > 0" :x="state.ins_xscale(0)" :y="state.ins_feat_lbl_y" fill="#ffffff">TSD</text>
      <line :x1="state.ins_xscale(state.polyx_x1)" :x2="state.ins_xscale(state.polyx_x2)" :y1="state.ins_feat_y" :y2="state.ins_feat_y" stroke-width="4" stroke="#a0ffff"/>
      <text :x="state.ins_xscale(state.polyx_x1)" :y="state.ins_feat_lbl_y" fill="#ffffff">{{ state.polyx_x2 == state.ins_len ? "polyA" : "polyT" }}</text>
      </svg>
  </div>
</template>

<style scoped>
.xaxis {
  font-size: 18px;
}

.ins_label {
  font-size: 18px;
  text-anchor: middle;
}

.me_label {
  font-size: 18px;
}

h1 {
  font-weight: 500;
  font-size: 2.6rem;
  top: -10px;
}

h3 {
  font-size: 1.2rem;
}

.greetings h1,
.greetings h3 {
  text-align: center;
}

@media (min-width: 1024px) {
  .greetings h1,
  .greetings h3 {
    text-align: left;
  }
}
</style>
