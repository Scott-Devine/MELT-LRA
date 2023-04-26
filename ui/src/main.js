import { createApp } from 'vue'

import App from './App.vue'

// Vuetify
import 'vuetify/styles'
import { createVuetify } from 'vuetify'
import * as components from 'vuetify/components'
import * as directives from 'vuetify/directives'

// D3
import { select } from 'd3-selection';
import 'd3-transition';

const vuetify = createVuetify({
  components,
  directives
})

import './assets/main.css'

const app = createApp(App).use(vuetify)
app.directive('axis', (el, binding) => {
    const axisMethod = binding.value
    select(el)
      .transition()
      .call(axisMethod)
})

app.mount('#app')
