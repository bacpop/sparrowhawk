<script setup lang="ts">
import type { ColumnDef, SortingState, ExpandedState } from '@tanstack/vue-table'
import {
    FlexRender,
    getCoreRowModel,
    getSortedRowModel,
    useVueTable,
} from '@tanstack/vue-table'
import { ref } from 'vue'
import { valueUpdater } from '@/lib/utils'
import type { TaxonomicIDRow } from './columns'
import { aniTextClass } from './columns'
import {
    Table,
    TableBody,
    TableCell,
    TableHead,
    TableHeader,
    TableRow,
} from '@/components/ui/table'

const props = defineProps<{
    columns: ColumnDef<TaxonomicIDRow>[]
    data: TaxonomicIDRow[]
}>()

const sorting = ref<SortingState>([])
const expanded = ref<ExpandedState>({})

const table = useVueTable({
    get data() { return props.data },
    get columns() { return props.columns },
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getRowCanExpand: row => (row.original.clusterDetails?.length ?? 0) > 0,
    onSortingChange: updaterOrValue => valueUpdater(updaterOrValue, sorting),
    onExpandedChange: updaterOrValue => valueUpdater(updaterOrValue, expanded),
    state: {
        get sorting() { return sorting.value },
        get expanded() { return expanded.value },
    },
})
</script>

<template>
    <div class="rounded-md border overflow-hidden">
        <Table>
            <TableHeader>
                <TableRow v-for="headerGroup in table.getHeaderGroups()" :key="headerGroup.id">
                    <TableHead v-for="header in headerGroup.headers" :key="header.id">
                        <FlexRender
                            v-if="!header.isPlaceholder"
                            :render="header.column.columnDef.header"
                            :props="header.getContext()"
                        />
                    </TableHead>
                </TableRow>
            </TableHeader>
            <TableBody>
                <template v-if="table.getRowModel().rows?.length">
                    <template v-for="row in table.getRowModel().rows" :key="row.id">
                        <TableRow
                            :class="[
                                row.getCanExpand() ? 'cursor-pointer hover:bg-muted/40' : '',
                                row.getIsExpanded() ? 'bg-muted/30' : '',
                            ]"
                            @click="row.getCanExpand() && row.toggleExpanded()"
                        >
                            <TableCell v-for="cell in row.getVisibleCells()" :key="cell.id">
                                <FlexRender :render="cell.column.columnDef.cell" :props="cell.getContext()" />
                            </TableCell>
                        </TableRow>
                        <TableRow v-if="row.getIsExpanded()">
                            <TableCell :colspan="columns.length" class="bg-muted/20 p-4">
                                <div class="w-fit max-w-full rounded-md border bg-white overflow-x-auto">
                                    <table class="w-auto text-sm">
                                        <thead class="bg-gray-50 text-gray-600">
                                            <tr>
                                                <th class="px-3 py-2 text-left font-medium whitespace-nowrap">Rank</th>
                                                <th class="px-3 py-2 text-left font-medium whitespace-nowrap">Species per cluster</th>
                                                <th class="px-3 py-2 text-right font-medium whitespace-nowrap">ANI</th>
                                                <th class="px-3 py-2 text-left font-medium whitespace-nowrap">Gemsparcl ID</th>
                                                <th class="px-3 py-2 text-left font-medium">GTDB species composition</th>
                                            </tr>
                                        </thead>
                                        <tbody>
                                            <tr
                                                v-for="candidate in row.original.clusterDetails"
                                                :key="row.id + '-' + candidate.rank"
                                                class="border-t"
                                            >
                                                <td class="px-3 py-2 whitespace-nowrap">{{ candidate.rank }}</td>
                                                <td class="px-3 py-2 italic whitespace-nowrap">{{ candidate.species }}</td>
                                                <td :class="['px-3 py-2 text-right font-medium whitespace-nowrap', aniTextClass(candidate.ani)]">{{ (candidate.ani * 100).toFixed(1) }}%</td>
                                                <td class="px-3 py-2 whitespace-nowrap">{{ candidate.metaGemsparcl }}</td>
                                                <td class="px-3 py-2 max-w-xl whitespace-normal">{{ candidate.metaGtdb }}</td>
                                            </tr>
                                        </tbody>
                                    </table>
                                </div>
                            </TableCell>
                        </TableRow>
                    </template>
                </template>
                <TableRow v-else>
                    <TableCell :colspan="columns.length" class="h-24 text-center">
                        No results.
                    </TableCell>
                </TableRow>
            </TableBody>
        </Table>
    </div>
</template>
