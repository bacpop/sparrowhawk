<template>
  <SidebarProvider class="h-svh">
    <Sidebar collapsible="none" class="mt-3.5 mb-6 mx-6 min-w-[185px] rounded-xl flex flex-col">
      <SidebarHeader>
        <img src="sparrowhawk-transparent.png" alt="Sparrowhawk" class="h-[120px] w-fit">
      </SidebarHeader>
      <SidebarContent>
        <SidebarGroup>
          <SidebarGroupContent>
            <SidebarMenu>
              <SidebarMenuItem v-for="item in tabs"
                               :key="item.sidebar_label"
                               :class="item.comingSoon
                                 ? 'py-2 px-3 cursor-not-allowed'
                                 : (item.id === tabName
                                     ? 'py-2 px-3 cursor-pointer bg-white rounded-sm shadow-sm'
                                     : 'py-2 px-3 cursor-pointer')">
                <SidebarMenuButton
                  @click="!item.comingSoon && changeTab(item.id)"
                  :class="item.comingSoon
                    ? 'p-0 pointer-events-none text-gray-400'
                    : 'p-0 hover:bg-transparent cursor-pointer'">
                  <component :is="item.icon"/>
                  <span class="text-md">{{ item.sidebar_label }}</span>
                </SidebarMenuButton>
              </SidebarMenuItem>
            </SidebarMenu>
          </SidebarGroupContent>
        </SidebarGroup>
      </SidebarContent>
      <SidebarFooter>
        <div class="flex flex-col opacity-80">
          <Button variant="link" size="sm" class="cursor-pointer items-center justify-start p-0"
                  @click="changeTab('faq')">
            FAQ
          </Button>
        </div>
      </SidebarFooter>
    </Sidebar>

    <!-- Main Content -->
    <main class="bg-white mt-6 mb-6 rounded-tl-xl rounded-bl-xl border border-gray-200 border-r-0 flex-1 overflow-y-auto p-8">
      <SidebarTrigger class="md:hidden fixed top-3 left-3 z-50 bg-white shadow-md rounded-md" />
      <div v-if="tabName === 'Assembly'">
        <AssemblyPage :tabName="tabName">
          <KmerHistogram class="mt-6"/>
        </AssemblyPage>
      </div>

      <div v-else-if="tabName === 'TaxonomicID'">
        <TaxonomicIDPage :tabName="tabName"/>
      </div>

      <div v-else-if="tabName === 'Mapping'">
        <MappingAlignmentPage :tabName="tabName">
<!--          <template #mapping>-->
<!--            <ResultsDisplayMapping/>-->
<!--          </template>-->
        </MappingAlignmentPage>
      </div>

      <div v-else-if="tabName === 'Alignment'">
        <MappingAlignmentPage :tabName="tabName">
          <template #alignment>
            <ResultsDisplayAlignment class="mt-6 mx-6 border-t border-t-gray-200 pt-6"/>
          </template>
        </MappingAlignmentPage>
      </div>

      <div v-else-if="tabName === 'GeneCalling'">
        <GeneCallingPage :tabName="tabName"/>
      </div>

      <div v-else-if="tabName === 'HostDepletion'">
        <HostDepletionPage :tabName="tabName"/>
      </div>

      <div v-else-if="tabName === 'faq'">
        <FaqPage/>
      </div>
    </main>
  </SidebarProvider>
</template>

<script lang="ts">
import {defineComponent} from 'vue';
import {useStore} from 'vuex';
// eslint-disable-next-line
import {Codesandbox, ScanFace, TextAlignCenter, TreePine, Dna, Funnel, Pill} from "lucide-vue-next";

import AssemblyPage from './components/pages/AssemblyPage.vue';
import MappingAlignmentPage from './components/pages/MappingAlignmentPage.vue';
import TaxonomicIDPage from './components/pages/TaxonomicIDPage.vue';
import GeneCallingPage from "./components/pages/GeneCallingPage.vue";
import HostDepletionPage from "./components/pages/HostDepletionPage.vue";
import ResultsDisplayMapping from './components/ResultsDisplayMapping.vue';
import ResultsDisplayAlignment from './components/ResultsDisplayAlignment.vue';
import KmerHistogram from './components/KmerHistogram.vue';
import WorkerAssembler from '@/workers/Assembler.worker';
import WorkerMapper from '@/workers/Mapper.worker';
import WorkerDepleter from '@/workers/Depleter.worker';
import "@fontsource/ibm-plex-sans";
import {
  Sidebar,
  SidebarContent,
  SidebarFooter,
  SidebarGroup,
  SidebarGroupContent,
  SidebarHeader,
  SidebarMenu,
  SidebarMenuButton,
  SidebarMenuItem,
  SidebarProvider,
  SidebarTrigger
} from "@/components/ui/sidebar";
import {Button} from "@/components/ui/button";
import FaqPage from "@/components/pages/FaqPage.vue";

interface Tab {
  id: string;
  sidebar_label: string;
  label: string;
  icon: string;
  comingSoon?: boolean;
}

export default defineComponent({
  name: 'App',

  components: {
    FaqPage,
    Sidebar,
    SidebarFooter,
    SidebarHeader,
    SidebarMenuItem,
    SidebarMenuButton,
    SidebarMenu,
    SidebarGroupContent,
    SidebarGroup,
    SidebarContent,
    SidebarProvider,
    SidebarTrigger,
    Codesandbox,
    Dna,
    Funnel,
    TextAlignCenter,
    TreePine,
    ScanFace,
    Pill,
    AssemblyPage,
    MappingAlignmentPage,
    TaxonomicIDPage,
    GeneCallingPage,
    HostDepletionPage,
    KmerHistogram,
    ResultsDisplayMapping,
    ResultsDisplayAlignment,
    Button
  },

  setup() {
    const store = useStore();
    return {store};
  },

  data() {
    return {
      tabName: 'Assembly' as string,
      tabs: [
        {id: 'Assembly', sidebar_label: "Assembly", label: 'Assembly', icon: 'Codesandbox'},
        {id: 'Mapping', sidebar_label: "Mapping", label: 'Mapping (within species)', icon: 'TextAlignCenter'},
        {id: 'Alignment', sidebar_label: "Alignment", label: 'Alignment (within species)', icon: 'TreePine'},
        {id: 'TaxonomicID', sidebar_label: "Taxonomic ID", label: 'Taxonomic ID', icon: 'ScanFace'},
        {id: 'GeneCalling', sidebar_label: "Gene calling", label: 'Gene calling', icon: 'Dna'},
        {id: 'HostDepletion', sidebar_label: "Host depletion", label: 'Host depletion', icon: 'Funnel'},
        {id: 'amrdetection', sidebar_label: "AMR detection (soon)", label: 'AMR detection (soon)', icon: 'Pill', comingSoon: true} as Tab,
      ] as Tab[]
    }
  },

  mounted(): void {
    // Set initial tab from URL hash
    const hash = window.location.hash.slice(1); // Remove the '#'
    const validTabs = this.tabs.filter(t => !t.comingSoon).map(t => t.id).concat(['faq']);
    if (hash && validTabs.includes(hash)) {
      this.tabName = hash;
    } else {
      // Set default hash if none or invalid
      window.location.hash = this.tabName;
    }

    // Listen for back/forward navigation
    window.addEventListener('hashchange', () => {
      const newHash = window.location.hash.slice(1);
      if (newHash && validTabs.includes(newHash)) {
        this.tabName = newHash;
      }
    });

    console.log("Loading wasm modules in workers...")

    import("@/pkg")
        .then(() => {
          if (window.Worker) {
            const worker = new WorkerAssembler();
            this.store.commit('SET_WORKER', worker);
          } else {
            throw new Error("WebWorkers are not supported by this web browser.");
          }
        });
    import("@/pkg_ska")
        .then(() => {
          if (window.Worker) {
            const worker = new WorkerMapper();
            this.store.commit('SET_WORKER_SKA', worker);
          } else {
            throw new Error("WebWorkers are not supported by this web browser.");
          }
        });
    import("@/pkg_orphos-bridge")
        .then(() => {
          if (window.Worker) {
            this.store.dispatch('initCallerWorkers', 4);
          } else {
            throw new Error("WebWorkers are not supported by this web browser.");
          }
        });
    import("@/pkg_sketchlib")
        .then(() => {
            if (window.Worker) {
                this.store.dispatch('initSketchlibWorkers', 4);
            } else {
                throw new Error("WebWorkers are not supported by this web browser.");
            }
       });
    import("@/pkg_deacon")
        .then(() => {
            if (window.Worker) {
                const worker = new WorkerDepleter();
                this.store.commit('SET_WORKER_DEACON', worker);
            } else {
                throw new Error("WebWorkers are not supported by this web browser.");
            }
        });
  },

  methods: {
    changeTab(tab: string): void {
      this.tabName = tab;
      window.location.hash = tab;
    }
  }
});

</script>
