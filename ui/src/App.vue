<script setup>
  import { ref, computed, reactive, watch } from 'vue';
  import MEI_Viewer from './components/MEI_Viewer.vue'
  import CALU_Report from './components/CALU_Report.vue'

  const csv_headers = ['chrom', 'pos', 'strand', 'ME', '%ME', '%id', '%id_ng', '%cov', 'insertion_seq', 'left_flank_seq', 'right_flank_seq', 'TSD_seq', 'polyX_coords', 'ME_coords', 'insertion_coords', 'match_string',
            'ME_family', 'ME_subfamily', 'ME_start', 'ME_stop', 'ME_num_diag_matches', 'ME_num_diffs', 'ME_diffs', 'overlapping_annots', 'genotype', 'hap1_region', 'hap2_region']

  function parse_mei(mei_str) {
    let f = mei_str.split(',')
    let nh = csv_headers.length
    let mei = {};
    for (let i = 0;i < nh; ++i) {
        if (f[i].endsWith('%')) {
            f[i] = f[i].slice(0, -1) * 1.0
        }
        mei[csv_headers[i]] = f[i];
    }
    // overlapping annotations/repeats
    let annots = mei['overlapping_annots'].split(/\|/).filter(a => a != '')
    mei['annots'] = annots
    mei['n_overlapping_annots'] = annots.length

    // unique overlapping repeat family/families
    let overlapping_rep_fams = {}
    annots.map(a => { 
        let fam = a.split(/:/)[6]
        overlapping_rep_fams[fam]++
    })
    mei['overlapping_rep_fam'] = Object.keys(overlapping_rep_fams).sort().join(",")
    if (mei['overlapping_rep_fam'] == "") mei['overlapping_rep_fam'] = 'None'

    mei['pos'] = mei['pos'] * 1.0
    mei['insertion_size'] = mei['insertion_seq'].length
    mei['TSD_length'] = mei['TSD_seq'].length

    let pxc = mei['polyX_coords'].split('-')
    mei['polyX_length'] = pxc[1] - pxc[0] >= 0 ? pxc[1] - pxc[0] : 0
    mei['key'] = mei['chrom'] + ':' + mei['pos']
    return mei
  }

  // hard-coded data sources
  const mei_files = [
    // HG00514 chrY only - DEBUG
    "assets/data/HG00514-chr22-PAV-MEs-v1.3.0-90-90-95bp.csv",
    // HG00514 all
    "assets/data/HG00514-PAV-MEs-v1.3.0-90-90-95bp.csv",
    // HG00514 all but chrY
//    "assets/data/HG00514-CCS-PAV-MEs-v1.2.0-90-90-95bp.csv",
    // 2023-08-23 5 new samples
    "assets/data/HG00171-PAV-MEs-v1.3.0-90-90-95bp.csv",
    "assets/data/HG00733-PAV-MEs-v1.3.0-90-90-95bp.csv",
    "assets/data/HG02666-PAV-MEs-v1.3.0-90-90-95bp.csv",
    "assets/data/HG02953-PAV-MEs-v1.3.0-90-90-95bp.csv",
    "assets/data/NA21487-PAV-MEs-v1.3.0-90-90-95bp.csv"
  ]

  const mei_urls = mei_files.map(f => { return new URL(f, import.meta.url).href })
    
  const state = reactive({
    'mei_urls': mei_urls,
    'selected_mei_url': null,
    'meis': [],
    'tab': "viewer"
  })

  watch(() => state.selected_mei_url, (newValue) => {
    fetch(newValue)
    .then(res => {
      return res.text()
    }).then(txt => {
      state.meis = txt.split("\n").slice(1).filter(l => !l.match(/^\s*$/)).map(ms => parse_mei(ms))
    })
    .catch(err => {
      console.log("caught error " + err)
    })
  })

  state.selected_mei_url = mei_urls[0]

</script>

<template>
  <v-app class="pa-0 ma-0 mr-3">
    <v-toolbar density="compact" color="primary" app>
      <v-app-bar-nav-icon></v-app-bar-nav-icon>
      <v-toolbar-title>
        MEI Callset: 
        <v-menu><template v-slot:activator="{ props }"><v-chip size="large" v-bind="props">{{state.selected_mei_url}}</v-chip></template><v-list v-model="state.mei_list"><v-list-item v-for="(i, ind) in state.mei_urls" @click="state.selected_mei_url = i"><v-list-item-title>{{i}}</v-list-item-title></v-list-item></v-list></v-menu>
      </v-toolbar-title>
      <template v-slot:extension>
        <v-tabs v-model="state.tab">
          <v-tab key="viewer" value="viewer">MEI Viewer</v-tab>
          <v-tab key="report" value="report">CALU/LINEU</v-tab>
        </v-tabs>
      </template>
    </v-toolbar>
    <v-main class="pa-0 ma-0" color="white" app>
      <v-container fluid class="pa-0 ma-0">
        
        <v-card-text>
          <v-window v-model="state.tab">
            <v-window-item key="viewer" value="viewer">
               <MEI_Viewer :meis="state.meis" />
            </v-window-item>
            <v-window-item key="report" value="report">
               <CALU_Report :meis="state.meis"/>
            </v-window-item>
          </v-window>
        </v-card-text>
      </v-container>
    </v-main>
  </v-app>
</template>


<style scoped>
.v-tab {
  font-size: 1.1rem;
  font-weight: bold !important;
}
</style>
