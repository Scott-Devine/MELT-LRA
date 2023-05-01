<script setup>

import { reactive } from 'vue'
import MEI from './MEI.vue'
import MiniMEI from './MiniMEI.vue'

const headers = [
    { text: 'chrom', value: 'chrom', sortable: true, fixed: true },
    { text: 'pos', value: 'pos', sortable: true, fixed: true },
    { text: 'strand', value: 'strand', sortable: true, fixed: true },
    { text: 'ME', value: 'ME', sortable: true, fixed: true },
    { text: 'size', value: 'insertion_size', sortable: true, fixed: true },
    { text: '%ME', value: '%ME', sortable: true, fixed: true },
    { text: '%id', value: '%id', sortable: true, fixed: true },
    { text: '%id_ng', value: '%id_ng', sortable: true, fixed: true },
    { text: '%cov', value: '%cov', sortable: true, fixed: true },
    { text: 'TSD_seq', value: 'TSD_seq', sortable: true, fixed: true},
    { text: 'TSD_bp', value: 'TSD_length', sortable: true, fixed: true},
    { text: 'polyX_bp', value: 'polyX_length', sortable: true, fixed: true},
    { text: 'MEI', value: 'MEI', fixed: true, width: 360 }
]

const state = reactive({
    meis: [],
    headers: headers,
    pctid: 90.0
})

const csv_headers = ['chrom', 'pos', 'strand', 'ME', '%ME', '%id', '%id_ng', '%cov', 'insertion_seq', 'left_flank_seq', 'right_flank_seq', 'TSD_seq', 'polyX_coords', 'ME_coords', 'insertion_coords', 'match_string']
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
//const mei_file = "../assets/data/chr22-MEs-90-90-95bp-v1.1.1.csv"
const mei_file = "../assets/data/HG00514-CCS-PAV-MEs-90-90-95bp-v1.1.1.csv"
const mei_url = new URL(mei_file, import.meta.url).href

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
    <v-app-bar color="info" :title="mei_url">
    </v-app-bar>

    <v-main class="px-0 mx-0">
        <v-container fluid style="width: 100vw;"></v-container>
        <v-card class="pa-4">
        
            <v-slider v-model="state.pctid" label="Min %identity:" :min="90" :max="100" :step="1" thumb-label>
            <template v-slot:append>
                <v-text-field v-model="state.pctid" hide-details type="number" style="width: 120px;"></v-text-field>
            </template>
        </v-slider>

        <v-chip size="x-large" color="primary">? / {{ state.meis.length }}</v-chip> <span class="text-subtitle-1">MEIs selected</span>
     
        <EasyDataTable
           :headers="state.headers"
            :items="state.meis"
            alternating
            border-cell
            buttons-pagination
            show-index
            style="width: 100%;"
            class="mt-4"
        >
            <template #item-MEI="item">
                <MiniMEI :key="item.key" :mei="item" />
            </template>

            <template #expand="item">
                <div class="px-2">
                    <MEI :key="item.key" :mei="item" />
                </div>
            </template>
        </EasyDataTable>

    </v-card>

<!--    <div v-for="(mei, m) in state.meis">
        <MEI v-if="mei['%id_ng'] >= state.pctid" :mei="mei" :label="(m+1) + ' / ' + state.meis.length" />
    </div>
-->
    </v-main>
</template>
  
<style scoped></style>