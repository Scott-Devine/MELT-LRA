<script setup>

import { ref, computed, reactive, watch } from 'vue';
import MEI from './MEI.vue'
import MiniMEI from './MiniMEI.vue'
import { format } from 'd3-format'

const props = defineProps({
    meis: Array
});

const pct_format = format(".1f")

const headers = [
{ text: 'chrom', value: 'chrom', sortable: true, fixed: true },
{ text: 'pos', value: 'pos', sortable: true, fixed: true },
{ text: 'strand', value: 'strand', sortable: true, fixed: true },
{ text: 'ME_family', value: 'ME_family', sortable: true, fixed: true },
{ text: 'ME_subfamily', value: 'ME_subfamily', sortable: true, fixed: true },
{ text: 'ins_size', value: 'insertion_size', sortable: true, fixed: true },
{ text: '%ME', value: '%ME', sortable: true, fixed: true },
{ text: '%id', value: '%id', sortable: true, fixed: true },
{ text: '%id_nogaps', value: '%id_ng', sortable: true, fixed: true },
{ text: '%coverage', value: '%cov', sortable: true, fixed: true },
{ text: 'TSD_bp', value: 'TSD_length', sortable: true, fixed: true},
{ text: 'polyA/T_bp', value: 'polyX_length', sortable: true, fixed: true},
{ text: 'MEI', value: 'MEI', fixed: true, width: 400 }
]

// EasyDataTable pagination
const dataTable = ref()
const nextPage = () => { dataTable.value.nextPage() }
const prevPage = () => { dataTable.value.prevPage() }
const gotoPage = (p) => { if (dataTable.value != null) dataTable.value.updatePage(p) }
const maxPage = computed(() => dataTable.value?.maxPaginationNumber);
const currentPage = computed(() => dataTable.value?.currentPaginationNumber);

const sortBy = []
const sortType = []
const me_types = ['ALU', 'LINE1', 'SVA']
const me_families = ['AluJ', 'AluS', 'AluY', 'LINE1', 'SVA']

const state = reactive({
    meis: [],
    total_mei_counts: {'ALU': 0, 'LINE1':0, 'SVA': 0},
    selected_meis: [],
    selected_mei_counts: {'ALU': 0, 'LINE1':0, 'SVA': 0},
    selected_me_types: ['ALU', 'LINE1', 'SVA'],
    selected_me_families: ['AluJ', 'AluS', 'AluY', 'LINE1', 'SVA'],
    headers: headers,
    pctid_range: [0.0, 100.0],
    min_pctid_nogaps: 90.0,
    min_pctcov: 90.0,
    me_pctcov_range: [0.0, 100.0],
    me_ins_length_range: [0.0, 7000.0],
    tsd_length_range: [0.0, 500.0],
    polyx_length_range: [0.0, 500.0],
    display_mode: 'table'
})

function update_selected_meis() {
    gotoPage(1)
    let f_meis = []

    // selected by type
    let selected_by_type = {}
    let n_selected_by_type = {}
    me_types.forEach(mt => {
        n_selected_by_type[mt] = 0
        selected_by_type[mt] = false
    })
    state.selected_me_types.forEach(mt =>{
        selected_by_type[mt] = true
    })

    // selected by family
    let selected_by_family = {}
    let n_selected_by_family = {}
    me_families.forEach(mf => {
        n_selected_by_family[mf] = 0
        selected_by_family[mf] = false
    })
    state.selected_me_families.forEach(mf =>{
        selected_by_family[mf] = true
    })

    state.meis.forEach(m => {
        if ((m['%id_ng'] >= state.min_pctid_nogaps) 
        && (m['%cov'] >= state.min_pctcov) 
        && (m['%id'] >= state.pctid_range[0]) && (m['%id'] <= state.pctid_range[1])
        && (m['%ME'] >= state.me_pctcov_range[0]) && (m['%ME'] <= state.me_pctcov_range[1])
        && (m['insertion_size'] >= state.me_ins_length_range[0]) && (m['insertion_size'] <= state.me_ins_length_range[1])
        && (m['TSD_length'] >= state.tsd_length_range[0]) && (m['TSD_length'] <= state.tsd_length_range[1])
        && (m['polyX_length'] >= state.polyx_length_range[0]) && (m['polyX_length'] <= state.polyx_length_range[1])
        && (selected_by_type[m['ME']])
        && (selected_by_family[m['ME_family']])
        ) {
            f_meis.push(m)
            n_selected_by_type[m.ME]++
            n_selected_by_family[m.ME_family]++
        }
    })
    state.selected_meis = f_meis
    state.selected_mei_counts = n_selected_by_type
    if (state.selected_meis.length == state.meis.length) state.total_mei_counts = n_selected_by_type
}

watch(() => props.meis, (newValue) => { state.meis = newValue })
watch(() => state.min_pctid_nogaps, (newValue) => { update_selected_meis() })
watch(() => state.pctid_range, (newValue) => { update_selected_meis() })
watch(() => state.min_pctcov, (newValue) => { update_selected_meis() })
watch(() => state.me_pctcov_range, (newValue) => { update_selected_meis() })
watch(() => state.me_ins_length_range, (newValue) => { update_selected_meis() })
watch(() => state.tsd_length_range, (newValue) => { update_selected_meis() })
watch(() => state.polyx_length_range, (newValue) => { update_selected_meis() })
watch(() => state.meis, (newValue) => { update_selected_meis() })
watch(() => state.selected_me_types, (newValue) => { update_selected_meis() })
watch(() => state.selected_me_families, (newValue) => { update_selected_meis() })

function formatRatio(n1,n2) {
    return String(n1).padStart(String(n2).length, " ") + "/" + n2 + " (" + pct_format((n1/n2) * 100.0) + "%)";
}

function getCountRatio(m) {
    const n1 = state.selected_mei_counts[m];
    const n2 = state.total_mei_counts[m];
    return formatRatio(n1, n2)
}

function getMEColor(me) {
    if (me == 'ALU') return "#1b9e77"
    if (me == 'LINE1') return "#d95f02"
    return "#7570b3";
}

state.meis = props.meis

</script>

<template>
    <v-card>
        <v-card variant="outlined" class="pa-0 ma-0">
            <v-card-title>
                <v-icon large class="pr-2">mdi-tune</v-icon><span class="font-weight-medium mr-4">Filter MEIs:</span>
                <v-chip label size="large" color="black" class="font-weight-medium ml-4">{{ formatRatio(state.selected_meis.length, state.meis.length) }}</v-chip> 
                <span class="text-h6 ml-2">total </span>
                <span v-for="me_type in ['ALU', 'LINE1', 'SVA']" class="text-h6 ml-3 mr-4">
                    <v-chip label size="large" :disabled="!state.selected_mei_counts[me_type]" :color="getMEColor(me_type)" class="font-weight-medium">{{ getCountRatio(me_type) }}</v-chip>
                    {{ me_type }}</span>
                    
                </v-card-title>
                <v-container class="pa-0 ma-0 pt-2 pl-4">
                    <v-row class="pa-0 ma-0">
                        <v-col cols="12" class="pa-0 ma-0">
                            <v-container class="pa-0 ma-0">
                                <v-row v-if="false" class="pa-0 ma-0">
                                    <v-col cols="3" class="pa-0 ma-0">
                                        ME type(s):
                                    </v-col>
                                    <v-col cols="6" class="pa-0 ma-0">
                                        <v-select v-model="state.selected_me_types" :items="['ALU', 'LINE1', 'SVA']" multiple hide-details variant="outlined" density="compact" class="pa-0 ma-0 pb-2"></v-select>
                                    </v-col>
                                    <v-col cols="3" class="pa-0 ma-0 pl-3">
                                        
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="3" class="pa-0 ma-0">
                                        ME families:
                                    </v-col>
                                    <v-col cols="6" class="pa-0 ma-0">
                                        <v-select v-model="state.selected_me_families" :items="me_families" multiple hide-details variant="outlined" density="compact" class="pa-0 ma-0 pb-2"></v-select>
                                    </v-col>
                                    <v-col cols="3" class="pa-0 ma-0 pl-3">
                                        
                                    </v-col>
                                </v-row>

                                <v-row class="pa-0 ma-0">
                                    <v-col cols="3" class="pa-0 ma-0">
                                        Insertion size range:
                                    </v-col>
                                    <v-col cols="6" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.me_ins_length_range" :min="0" :max="7000" :step="50" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="3" class="pa-0 ma-0 pl-3">
                                        {{ state.me_ins_length_range[0] }}bp - {{ state.me_ins_length_range[1] }}bp
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="3" class="pa-0 ma-0">
                                        Percent identity range:
                                    </v-col>
                                    <v-col cols="6" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.pctid_range" :min="0" :max="100" :step="1" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="3" class="pa-0 ma-0 pl-3">
                                        {{ state.pctid_range[0] }}% - {{ state.pctid_range[1] }}%
                                    </v-col>
                                </v-row>

                                <v-row class="pa-0 ma-0">
                                    <v-col cols="3" class="pa-0 ma-0">
                                        Min ungapped %identity:
                                    </v-col>
                                    <v-col cols="6" class="pa-0 ma-0">
                                        <v-slider v-model="state.min_pctid_nogaps" :min="90" :max="100" :step="1" thumb-label hide-details></v-slider>
                                    </v-col>
                                    <v-col cols="3" class="pa-0 ma-0 pl-3">
                                        {{ state.min_pctid_nogaps }}%
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="3" class="pa-0 ma-0">
                                        Min insertion %coverage:
                                    </v-col>
                                    <v-col cols="6" class="pa-0 ma-0">
                                        <v-slider v-model="state.min_pctcov" :min="90" :max="100" :step="1" thumb-label hide-details></v-slider>
                                    </v-col>
                                    <v-col cols="3" class="pa-0 ma-0 pl-3">
                                        {{ state.min_pctcov }}%
                                    </v-col>
                                </v-row>

                                <v-row class="pa-0 ma-0">
                                    <v-col cols="3" class="pa-0 ma-0">
                                        Reference ME %coverage range:
                                    </v-col>
                                    <v-col cols="6" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.me_pctcov_range" :min="0" :max="100" :step="1" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="3" class="pa-0 ma-0 pl-3">
                                        {{ state.me_pctcov_range[0] }}% - {{ state.me_pctcov_range[1] }}%
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="3" class="pa-0 ma-0">
                                        TSD length range:
                                    </v-col>
                                    <v-col cols="6" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.tsd_length_range" :min="0" :max="500" :step="1" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="3" class="pa-0 ma-0 pl-3">
                                        {{ state.tsd_length_range[0] }}bp - {{ state.tsd_length_range[1] }}bp
                                    </v-col>
                                </v-row>
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="3" class="pa-0 ma-0">
                                        polyA/T length range:
                                    </v-col>
                                    <v-col cols="6" class="pa-0 ma-0">
                                        <v-range-slider v-model="state.polyx_length_range" :min="0" :max="500" :step="1" thumb-label hide-details></v-range-slider>
                                    </v-col>
                                    <v-col cols="3" class="pa-0 ma-0 pl-3">
                                        {{ state.polyx_length_range[0] }}bp - {{ state.polyx_length_range[1] }}bp
                                    </v-col>
                                </v-row>         
                                
                                <v-row class="pa-0 ma-0">
                                    <v-col cols="3" class="pa-0 ma-0">
                                        Display:
                                    </v-col>
                                    <v-col cols="9" class="pa-0 ma-0">
                                        <v-radio-group v-model="state.display_mode" inline>
                                            <v-radio label="Sortable table" value="table" density="compact"></v-radio>
                                            <v-radio label="Full-size figures [much slower]" value="figures" density="compact" class="pl-2"></v-radio>
                                        </v-radio-group>
                                    </v-col>
                                </v-row>       
                            </v-container>
                            
                        </v-col>
                    </v-row>
                </v-container>
            </v-card>
            
            <!-- sortable table view -->
            <v-card v-if="state.display_mode == 'table'" class="pa-2">
                <!-- supplemental pagination controls -->
                <div class="pa-2" style="background-color: #e0e0e0;">
                    <v-btn @click="gotoPage(1)" density="compact" :disabled="maxPage == 0">Page 1</v-btn>
                    <v-btn @click="prevPage()" density="compact" prepend-icon="mdi-arrow-left-bold" class="ml-2" :disabled="maxPage == 0">previous page</v-btn>
                    <v-btn @click="nextPage()" density="compact" append-icon="mdi-arrow-right-bold" class="ml-2" :disabled="maxPage == 0">next page</v-btn>
                    <span v-if="maxPage > 0" class="px-3 font-weight-medium">Page {{ currentPage }} / {{ maxPage }}</span>
                    <span v-else>No MEIs selected</span>
                </div>
                <EasyDataTable
                ref="dataTable"
                :headers="state.headers"
                :items="state.selected_meis"
                alternating
                border-cell
                :sort-by="sortBy"
                :sort-type="sortType"
                multi-sort
                buttons-pagination
                show-index
                class="mt-1"
                >
                <template #item-MEI="item">
                    <MiniMEI :key="item.key + '-mini'" :mei="item" />
                </template>
                
                <template #expand="item">
                    <div class="px-2">
                        <MEI :key="item.key" :mei="item" />
                        <div class="text-h6">
                            <div class="tsd_div pa-1 my-1 font-weight-bold">TSD</div> {{ item.TSD_seq }}<br>
                            <span v-if="item.ME == 'ALU' || item.ME == 'LINE1'">
                                <div class="calu_div pa-1 my-1 font-weight-bold">{{ item.ME == 'ALU' ? 'CALU' : 'LINEU'}}</div> 
                                {{item.ME_family}} {{item.ME_subfamily}} {{item.ME_start}}-{{item.ME_stop}} diag_matches={{item.ME_num_diag_matches}} num_diffs={{ item.ME_num_diffs }}  diffs={{ item.ME_diffs }}<br>
                            </span>
                        </div>
                    </div>
                </template>           
            </EasyDataTable>
            <!-- supplemental pagination controls -->
            <div class="pa-2" style="background-color: #e0e0e0;">
                <v-btn @click="gotoPage(1)" density="compact" :disabled="maxPage == 0">Page 1</v-btn>
                <v-btn @click="prevPage()" density="compact" prepend-icon="mdi-arrow-left-bold" class="ml-2" :disabled="maxPage == 0">previous page</v-btn>
                <v-btn @click="nextPage()" density="compact" append-icon="mdi-arrow-right-bold" class="ml-2" :disabled="maxPage == 0">next page</v-btn>
                <span v-if="maxPage > 0" class="px-3 font-weight-medium">Page {{ currentPage }} / {{ maxPage }}</span>
                <span v-else>No MEIs selected</span>
            </div>
        </v-card>
        <!-- list of figures view -->
        <v-card v-else class="pa-2" style="background-color: black;">
            <div v-for="(mei, m) in state.selected_meis">
                <MEI :mei="mei" :key="mei.key" :label="(m + 1) + '/' + state.selected_meis.length" />
            </div>
        </v-card>
    </v-card>
</template>

<style scoped>
div.tsd_div {
    display: inline-block;
    background-color: #a0ffa0;
    border: 1px solid black;
}
div.calu_div {
    display: inline-block;
    background-color: #d0d0d0;
    border: 1px solid black;
}
</style>
