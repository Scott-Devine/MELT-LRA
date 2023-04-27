<script setup>

  import { reactive } from 'vue'
  import { format } from 'd3-format'
  import { scaleLinear, scaleLog } from 'd3-scale';
  import { axisTop, axisBottom, axisLeft, axisRight } from 'd3-axis';
  
  const props = defineProps({
    meiStr: { type: String, required: true }
  })
console.log("props=" + props)
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
  const mei = parse_mei(props.meiStr)
  
  // SVG coordinate system
  const width = 1200
  const height = 250
  const margins = { left: 30, right: 30, top: 40, bottom: 20}
  let rx1 = margins.left
  let rx2 = width - margins.right
 
  const left_flank_bp = mei.left_flank_seq.length
  const right_flank_bp = mei.right_flank_seq.length
  const flank_bp = left_flank_bp + right_flank_bp

  const ins_len = mei.insertion_seq.length
  const me_len = ME_LENGTHS[mei.ME]
  // center the two sequences horizontally
  const max_len = (ins_len > me_len ? ins_len : me_len) + flank_bp;
  const max_xscale = scaleLinear().domain([0, max_len]).range([rx1, rx2])
  console.log("ins_len=" + ins_len +  " me_len=" + me_len)

  // insertion coordinate system
  const ins_x = (max_len - ins_len)/2
  const ins_rx1 = max_xscale(ins_x)
  const ins_rx2 = max_xscale(ins_x + ins_len)
  const ins_xscale = scaleLinear().domain([0, ins_len]).range([ins_rx1, ins_rx2]);
  const ins_xaxis = axisTop(ins_xscale)
  const[px_x1, px_x2] = mei.polyX_coords.split("-").map(x => x * 1.0)
  const[ins_x1, ins_x2] = mei.insertion_coords.split("-").map(x => x * 1.0)

  // reference ME coordinate system
  const me_x = (max_len - me_len)/2
  const me_rx1 = max_xscale(me_x)
  const me_rx2 = max_xscale(me_x + me_len)
  const me_xscale = scaleLinear().domain([0, me_len]).range([me_rx1, me_rx2]);
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
    if (((ins_o2 - ins_o1) > 0) && ((me_o2 - me_o1) > 0)) {
      // TODO - handle reverse strand matches
      const span = {
        'ins': [ins_x1 + ins_o1 - 1, ins_x1 + ins_o2 - 1], 
        'me': [me_x1 + me_o1 - 1, me_x1 + me_o2 - 1],
        'pct_id': (n_id_bp / (ins_o2 - ins_o1)) * 100.0
      }
      spans.push(span)
      console.log("adding " + span.ins[0] + ' - ' + span.ins[1] + " / " + span.me[0] + ' - ' + span.me[1] + " " + span.pct_id)
      n_id_bp = 0
    }
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

  // convert spans to polygons
  spans.forEach(s => {
    s.points = []
    s.points.push([ins_xscale(s.ins[0]), 50])
    s.points.push([ins_xscale(s.ins[1]), 50])
    s.points.push([me_xscale(s.me[1]), 180])
    s.points.push([me_xscale(s.me[0]), 180])
    s.points_str = s.points.join(" ")
  })

  const state = reactive({
     // SVG
     width: width,
    height: height,
    margins: margins,

    mei: mei, 
    ins_len: ins_len,
    tsd_len: mei.TSD_seq.length,
    // insertion and associated features
    ins_xscale: ins_xscale, 
    ins_xaxis: ins_xaxis,
    polyx_x1: px_x1,
    polyx_x2: px_x2,
    ins_x1: ins_x1,
    ins_x2: ins_x2,
    // reference ME
    me_xscale: me_xscale,
    me_xaxis: me_xaxis,
    me_x1: me_x1,
    me_x2: me_x2,
    // alignment
    spans: spans
  })
  
</script>

<template>
  <div class="greetings">
    <h1 class="green">[{{state.mei.ME}}] {{ state.mei.chrom }}:{{ state.mei.pos }}:{{ state.mei.strand }}:{{state.ins_len}} bp</h1>

    <svg
      ref="mei_svg"
      key="mei-svg-1"
      :width="state.width"
      :height="state.height"
      class="pa-0 ma-0"
      style="border: 1px solid white"
      >

      <!-- insertion -->
      <g v-axis="state.ins_xaxis" class="xaxis" :transform="`translate(0,${state.margins.top})`">
      </g>
      <line v-if="state.tsd_len > 0" :x1="state.ins_xscale(0)" :x2="state.ins_xscale(state.tsd_len)" :y1="50" :y2="50" stroke-width="4" stroke="#a0ffa0"/>
      <text v-if="state.tsd_len > 0" :x="state.ins_xscale(0)" :y="70" fill="#ffffff">TSD</text>
      <line :x1="state.ins_xscale(state.polyx_x1)" :x2="state.ins_xscale(state.polyx_x2)" :y1="50" :y2="50" stroke-width="4" stroke="#a0ffff"/>
      <text :x="state.ins_xscale(state.polyx_x1)" :y="70" fill="#ffffff">polyX</text>
      <!-- <line :x1="state.ins_xscale(state.ins_x1)" :x2="state.ins_xscale(state.ins_x2)" :y1="100" :y2="10" stroke-width="4" stroke="#ffa0a0"/> -->

      <!-- reference ME -->
      <g v-axis="state.me_xaxis" class="xaxis" :transform="`translate(0,${state.height - 60})`">
      </g>

      <!-- matching spans -->
      <polygon v-for="s in state.spans" :points="s.points_str" fill="rgba(255,0,0,0.1)" stroke="#ff0000" />
      <text v-for="s in state.spans" :x="state.me_xscale((s.me[0] + s.me[1]) / 2) + 3" :y="160" fill="#ffffff" text-anchor="middle">{{ format('.0f')(s.pct_id) }}%</text>
      </svg>
  </div>
</template>

<style scoped>
.xaxis {
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
