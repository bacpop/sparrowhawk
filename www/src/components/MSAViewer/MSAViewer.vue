<script setup lang="ts">
import {computed, onMounted, onUnmounted, ref} from 'vue';
import {X} from "@lucide/vue";

// Simple throttle implementation
function throttle<T extends (...args: unknown[]) => void>(fn: T, wait: number): T {
  let lastTime = 0;
  let timeoutId: ReturnType<typeof setTimeout> | null = null;

  return function (this: unknown, ...args: Parameters<T>) {
    const now = Date.now();
    const remaining = wait - (now - lastTime);

    if (remaining <= 0) {
      if (timeoutId) {
        clearTimeout(timeoutId);
        timeoutId = null;
      }
      lastTime = now;
      fn.apply(this, args);
    } else if (!timeoutId) {
      timeoutId = setTimeout(() => {
        lastTime = Date.now();
        timeoutId = null;
        fn.apply(this, args);
      }, remaining);
    }
  } as T;
}

interface MSASequence {
  id: string;
  sequence: string;
}

interface Props {
  data: MSASequence[];
  filename?: string;
  showHeader?: boolean;
  initialLabelWidth?: number;
}

const props = withDefaults(defineProps<Props>(), {
  filename: 'Alignment',
  showHeader: true,
  initialLabelWidth: 100,
});

const emit = defineEmits<{
  remove: [];
}>();

// Constants
const CELL_SIZE = 24;
const RULER_HEIGHT = 16;
const OVERSCAN = 5;

// Color schemes - defined outside reactive system
const ResidueColorSchemes = {
  protein: {
    A: 'bg-green-200', C: 'bg-yellow-200', D: 'bg-red-200', E: 'bg-red-200',
    F: 'bg-purple-200', G: 'bg-gray-200', H: 'bg-pink-200', I: 'bg-blue-200',
    K: 'bg-orange-200', L: 'bg-blue-200', M: 'bg-blue-100', N: 'bg-red-100',
    P: 'bg-teal-200', Q: 'bg-red-100', R: 'bg-orange-300', S: 'bg-green-100',
    T: 'bg-green-100', V: 'bg-blue-100', W: 'bg-purple-300', Y: 'bg-purple-100',
    '-': 'bg-white',
  } as Record<string, string>,
  nucleotide: {
    'A': 'bg-green-200',
    'C': 'bg-blue-200',
    'G': 'bg-yellow-200',
    'T': 'bg-red-200',
    'U': 'bg-red-200',
    'N': 'bg-gray-200',
    '-': 'bg-white',
    '.': 'bg-white',
  } as Record<string, string>,
} as const;

const proteinOnlyChars = new Set(['D', 'E', 'F', 'H', 'I', 'K', 'L', 'M', 'P', 'Q', 'R', 'S', 'V', 'W', 'Y']);

// Refs
const scrollContainerRef = ref<HTMLElement | null>(null);
const labelWidth = ref(props.initialLabelWidth);
const scrollLeft = ref(0);
const scrollTop = ref(0);
const hoveredCol = ref<number | null>(null);
const hoveredRow = ref<number | null>(null);
const tooltipPos = ref({x: 0, y: 0});
const viewportWidth = ref(800);
const viewportHeight = ref(400);

// Drag state
const isDragging = ref(false);
const dragStartX = ref(0);

// Check if we have valid data
const hasData = computed(() => {
  const valid = props.data && props.data.length > 0 && props.data[0]?.sequence?.length > 0;
  console.log('[MSAViewer] hasData check:', {
    valid,
    dataLength: props.data?.length,
    firstSeqLength: props.data?.[0]?.sequence?.length,
  });
  return valid;
});

// Detect nucleotide vs protein - computed once per data change
const isNucleotide = computed(() => {
  if (!hasData.value) return true;
  for (let i = 0; i < Math.min(props.data.length, 10); i++) {
    const seq = props.data[i]?.sequence;
    if (!seq) continue;
    const upper = seq.toUpperCase();
    for (let j = 0; j < Math.min(upper.length, 50); j++) {
      if (proteinOnlyChars.has(upper[j])) return false;
    }
  }
  return true;
});

const colorScheme = computed(() =>
    isNucleotide.value ? ResidueColorSchemes.nucleotide : ResidueColorSchemes.protein
);

// Dimensions - simple computed values
const rowCount = computed(() => props.data?.length || 0);
const rawColCount = computed(() => props.data?.[0]?.sequence?.length || 0);
const colCount = computed(() => rawColCount.value);
const totalGridWidth = computed(() => colCount.value * CELL_SIZE);
const totalGridHeight = computed(() => rowCount.value * CELL_SIZE);

// Virtualization boundaries
const firstRow = computed(() => {
  if (!hasData.value) return 0;
  return Math.max(0, Math.floor(scrollTop.value / CELL_SIZE) - OVERSCAN);
});

const lastRow = computed(() => {
  if (!hasData.value) return 0;
  const visibleRows = Math.ceil(viewportHeight.value / CELL_SIZE);
  return Math.min(rowCount.value, firstRow.value + visibleRows + OVERSCAN * 2);
});

const firstCol = computed(() => {
  if (!hasData.value) return 0;
  return Math.max(0, Math.floor(scrollLeft.value / CELL_SIZE) - OVERSCAN);
});

const lastCol = computed(() => {
  if (!hasData.value) return 0;
  const visibleCols = Math.ceil((viewportWidth.value - labelWidth.value) / CELL_SIZE);
  return Math.min(colCount.value, firstCol.value + visibleCols + OVERSCAN * 2);
});

// Pre-computed cell data - this is the key optimization!
// We compute an array of cell objects with all their display properties
interface CellData {
  key: string;
  rowIdx: number;
  colIdx: number;
  char: string;
  bgClass: string;
  top: number;
  left: number;
}

const cellsData = computed((): CellData[] => {
  if (!hasData.value) {
    console.log('[MSAViewer] cellsData: no data');
    return [];
  }

  console.log('[MSAViewer] cellsData computing:', {
    firstRow: firstRow.value,
    lastRow: lastRow.value,
    firstCol: firstCol.value,
    lastCol: lastCol.value,
    rowCount: rowCount.value,
    colCount: colCount.value,
    rawColCount: rawColCount.value,
  });

  const cells: CellData[] = [];
  const scheme = colorScheme.value;
  const data = props.data;

  for (let r = firstRow.value; r < lastRow.value; r++) {
    const seq = data[r]?.sequence;
    if (!seq) continue;

    for (let c = firstCol.value; c < lastCol.value; c++) {
      const char = seq[c] || '';
      const bgClass = scheme[char.toUpperCase()] || 'bg-white';

      cells.push({
        key: `${r}-${c}`,
        rowIdx: r,
        colIdx: c,
        char,
        bgClass,
        top: r * CELL_SIZE,
        left: c * CELL_SIZE,
      });
    }
  }

  console.log('[MSAViewer] cellsData result:', cells.length, 'cells');
  return cells;
});

// Row labels data
interface RowLabelData {
  idx: number;
  id: string;
  top: number;
}

const rowLabelsData = computed((): RowLabelData[] => {
  if (!hasData.value) return [];

  const labels: RowLabelData[] = [];
  for (let r = firstRow.value; r < lastRow.value; r++) {
    labels.push({
      idx: r,
      id: props.data[r]?.id || '',
      top: r * CELL_SIZE,
    });
  }
  return labels;
});

// Ruler ticks data
interface RulerTickData {
  idx: number;
  left: number;
  label: string | null;
}

const rulerTicksData = computed((): RulerTickData[] => {
  if (!hasData.value) return [];

  const ticks: RulerTickData[] = [];
  const total = colCount.value;

  for (let c = firstCol.value; c < lastCol.value; c++) {
    const pos = c + 1;
    let label: string | null = null;

    if (pos % 10 === 0 || pos === 1 || pos === total) {
      label = String(pos);
    } else if (pos % 5 === 0) {
      label = '·';
    }

    if (label !== null) {
      ticks.push({
        idx: c,
        left: c * CELL_SIZE,
        label,
      });
    }
  }
  return ticks;
});

// Throttled scroll handler - critical for performance
const onScroll = throttle((e: Event) => {
  const target = e.target as HTMLElement;
  if (!target) return;
  scrollLeft.value = target.scrollLeft;
  scrollTop.value = target.scrollTop;
}, 16); // ~60fps

// Mouse move handler
function onMouseMove(e: MouseEvent) {
  if (!scrollContainerRef.value || !hasData.value) return;

  const rect = scrollContainerRef.value.getBoundingClientRect();
  const mouseX = e.clientX - rect.left;
  const mouseY = e.clientY - rect.top;

  const isOverCells = mouseX > labelWidth.value;

  const col = isOverCells
      ? Math.floor((mouseX - labelWidth.value + scrollLeft.value) / CELL_SIZE)
      : null;
  const row = Math.floor((mouseY - RULER_HEIGHT + scrollTop.value) / CELL_SIZE);

  hoveredRow.value = (row >= 0 && row < rowCount.value) ? row : null;
  hoveredCol.value = (col !== null && col >= 0 && col < colCount.value) ? col : null;
  tooltipPos.value = {x: e.clientX, y: e.clientY};
}

function onMouseLeave() {
  hoveredCol.value = null;
  hoveredRow.value = null;
}

// Label column drag resize
function onDragStart(e: MouseEvent) {
  isDragging.value = true;
  dragStartX.value = e.clientX;
  document.body.style.cursor = 'col-resize';
  document.addEventListener('mousemove', onDrag);
  document.addEventListener('mouseup', onDragEnd);
}

function onDrag(e: MouseEvent) {
  if (!isDragging.value) return;
  const delta = e.clientX - dragStartX.value;
  dragStartX.value = e.clientX;
  labelWidth.value = Math.max(40, Math.min(400, labelWidth.value + delta));
}

function onDragEnd() {
  isDragging.value = false;
  document.body.style.cursor = '';
  document.removeEventListener('mousemove', onDrag);
  document.removeEventListener('mouseup', onDragEnd);
}

// Resize observer
let resizeObserver: ResizeObserver | null = null;

onMounted(() => {
  if (scrollContainerRef.value) {
    resizeObserver = new ResizeObserver((entries) => {
      const entry = entries[0];
      if (entry) {
        viewportWidth.value = entry.contentRect.width;
        viewportHeight.value = entry.contentRect.height;
      }
    });
    resizeObserver.observe(scrollContainerRef.value);
  }
});

onUnmounted(() => {
  resizeObserver?.disconnect();
  document.removeEventListener('mousemove', onDrag);
  document.removeEventListener('mouseup', onDragEnd);
});
</script>

<template>
  <div class="relative flex flex-col w-full bg-white" style="max-width: 100%; contain: strict;">
    <!-- Header -->
    <div v-if="showHeader"
         class="bg-gradient-to-b from-gray-100 to-white px-2 py-1 flex items-center justify-between border-b">
      <span class="font-bold text-sm truncate flex-1 text-center">{{ filename }}</span>
      <button
          @click="emit('remove')"
          class="p-0.5 inline-flex items-center justify-center w-6 h-6 rounded-lg bg-gray-100 border border-gray-400 hover:bg-red-300"
      >
        <X class="w-4 h-4 text-red-700"/>
      </button>
    </div>

    <div v-if="!hasData" class="flex-1 flex items-center justify-center text-gray-400">
      No alignment data to display
    </div>


    <Teleport to="body">
      <div
          v-if="hoveredCol !== null && hasData"
          class="fixed px-2 py-1 text-sm bg-gray-100 rounded-xl pointer-events-none z-[9999] shadow border border-gray-400 font-sans"
          :style="{
            left: `${tooltipPos.x + 12}px`,
            top: `${tooltipPos.y + 12}px`,
          }"
      >
        <div class="flex flex-col items-center">
          <span class="font-bold">Site {{ hoveredCol + 1 }}</span>
          <span v-if="hoveredRow !== null" class="text-gray-600 text-xs">{{ data[hoveredRow]?.id }}</span>
        </div>
      </div>
    </Teleport>

    <div v-if="hasData" class="flex-1 relative min-h-0">
      <div
          ref="scrollContainerRef"
          class="absolute inset-0 overflow-auto font-mono text-sm"
          @scroll="onScroll"
          @mousemove="onMouseMove"
          @mouseleave="onMouseLeave"
      >
        <div
            class="relative"
            :style="{
              width: `${totalGridWidth + labelWidth}px`,
              height: `${totalGridHeight + RULER_HEIGHT}px`
            }"
        >
        <div class="sticky top-0 z-10 flex" :style="{ height: `${RULER_HEIGHT}px` }">
          <div
              class="sticky left-0 bg-white z-20"
              :style="{ width: `${labelWidth}px`, height: `${RULER_HEIGHT}px` }"
          />
          <div class="relative bg-white" :style="{ width: `${totalGridWidth}px` }">
            <div
                v-for="tick in rulerTicksData"
                :key="tick.idx"
                class="absolute flex items-center justify-center text-xs text-gray-600"
                :style="{
                  left: `${tick.left}px`,
                  width: `${CELL_SIZE}px`,
                  height: `${RULER_HEIGHT}px`,
                }"
            >
              {{ tick.label }}
            </div>
          </div>
        </div>

        <div
            class="sticky left-0 z-10 bg-white"
            :style="{
              position: 'absolute',
              top: `${RULER_HEIGHT}px`,
              width: `${labelWidth}px`,
              height: `${totalGridHeight}px`
            }"
        >
          <div
              v-for="label in rowLabelsData"
              :key="label.idx"
              :title="label.id"
              class="absolute flex items-center text-right font-bold truncate pl-2 pr-2 border-r border-gray-100"
              :class="{ 'bg-yellow-100': hoveredRow === label.idx }"
              :style="{
                top: `${label.top}px`,
                width: '100%',
                height: `${CELL_SIZE}px`,
                lineHeight: `${CELL_SIZE}px`,
              }"
          >
            <span class="block truncate">{{ label.id }}</span>
          </div>
          <div
              class="absolute top-0 right-0 w-2 h-full cursor-col-resize bg-black/5 hover:bg-black/10"
              @mousedown="onDragStart"
          />
        </div>

        <div
            :style="{
              position: 'absolute',
              top: `${RULER_HEIGHT}px`,
              left: `${labelWidth}px`,
              height: `${totalGridHeight}px`,
              width: `${totalGridWidth}px`
            }"
        >
          <div
              v-for="cell in cellsData"
              :key="cell.key"
              class="absolute flex items-center justify-center"
              :class="[cell.bgClass, { 'ring-1 ring-inset ring-blue-400': hoveredCol === cell.colIdx }]"
              :style="{
                top: `${cell.top}px`,
                left: `${cell.left}px`,
                height: `${CELL_SIZE}px`,
                width: `${CELL_SIZE}px`,
              }"
          >
            {{ cell.char }}
          </div>
        </div>
      </div>
    </div>
    </div>
  </div>
</template>
