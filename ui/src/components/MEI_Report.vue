<script setup>

import { ref, computed, reactive, watch } from 'vue';
import { findOrfs } from '../orfs.js'
import { symbolAsterisk } from 'd3';

const props = defineProps({
    meis: Array
});

const state = reactive({
    meis: []
})

function update_meis() {
    console.log("update_meis called, state has " + state.meis.length + " MEIs")
}

watch(() => props.meis, (newValue) => { console.log("props.meis updated"); state.meis = newValue })
watch(() => state.meis, (newValue) => { console.log("state.meis updated"); update_meis() })

state.meis = props.meis

</script>

<template>
    <v-container fluid class="mx-auto">
        <v-card variant="outlined" class="pa-2 ma-2 mx-auto">
            <v-card-title>
             title
                </v-card-title>
                <v-container class="pa-0 ma-0 pt-2 pl-4">
                    <v-row class="pa-0 ma-0">
                        <v-col cols="12" class="pa-0 ma-0">
                           
                        content goes here

                        </v-col>
                    </v-row>
                </v-container>
            </v-card>
        </v-container>
</template> 

<style scoped>
</style>
