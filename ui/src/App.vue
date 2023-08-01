<script setup>
  import { ref, computed, reactive, watch } from 'vue';
  import MEI_Viewer from './components/MEI_Viewer.vue'
  import MEI_Report from './components/MEI_Report.vue'

  const csv_headers = ['chrom', 'pos', 'strand', 'ME', '%ME', '%id', '%id_ng', '%cov', 'insertion_seq', 'left_flank_seq', 'right_flank_seq', 'TSD_seq', 'polyX_coords', 'ME_coords', 'insertion_coords', 'match_string',
            'ME_family', 'ME_subfamily', 'ME_start', 'ME_stop', 'ME_num_diag_matches', 'ME_num_diffs', 'ME_diffs']

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
    let pxc = mei['polyX_coords'].split('-')
    
    mei['pos'] = mei['pos'] * 1.0
    mei['insertion_size'] = mei['insertion_seq'].length
    mei['TSD_length'] = mei['TSD_seq'].length
    mei['polyX_length'] = pxc[1] - pxc[0] >= 0 ? pxc[1] - pxc[0] : 0
    mei['key'] = mei['chrom'] + ':' + mei['pos']
    return mei
  }

  // hard-coded data source
  //const mei_file = "../assets/data/chr22-MEs-v1.2.0-90-90-95bp.csv"
  const mei_file = "assets/data/HG00514-CCS-PAV-MEs-v1.2.0-90-90-95bp.csv"
  const mei_url = new URL(mei_file, import.meta.url).href

  const state = reactive({
    'meis': [],
    'tab': "viewer"
  })

  fetch(mei_url)
  .then(res => {
    return res.text()
  }).then(txt => {
    state.meis = txt.split("\n").slice(1).filter(l => !l.match(/^\s*$/)).map(ms => parse_mei(ms))
  })
  .catch(err => {
    console.log("caught error " + err)
  })

</script>

<template>
  <v-app class="pa-0 ma-0 mr-3">
    <v-toolbar density="compact" color="primary" app>
      <v-app-bar-nav-icon></v-app-bar-nav-icon>
      <v-toolbar-title>MEI callset: {{ mei_url }}</v-toolbar-title>
      <template v-slot:extension>
        <v-tabs v-model="state.tab">
          <v-tab key="viewer" value="viewer">MEI Viewer</v-tab>
          <v-tab key="report" value="report">MEI Summary</v-tab>
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
               <MEI_Report :meis="state.meis"/>
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
